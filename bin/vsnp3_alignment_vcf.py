#!/usr/bin/env python

__version__ = "3.16"

import os
import subprocess
import re
import glob
import shutil
import argparse
import textwrap
from datetime import datetime
import pandas as pd
from Bio import SeqIO

from vsnp3_file_setup import Setup
from vsnp3_file_setup import Banner
from vsnp3_file_setup import Latex_Report
from vsnp3_file_setup import Excel_Stats
from vsnp3_fastq_stats_seqkit import FASTQ_Stats
from vsnp3_vcf_annotation import VCF_Annotation
from vsnp3_assembly import Assemble
from vsnp3_zero_coverage import Zero_Coverage


class Alignment(Setup):
    ''' 
    '''

    def __init__(self, SAMPLE_NAME=None, FASTQ_R1=None, FASTQ_R2=None, reference=None, nanopore=False, gbk=None, assemble_unmap=None, debug=False):
        '''
        Start at class call
        '''
        Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=reference, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.print_run_time('Align and make VCF file')
        self.nanopore = nanopore
        self.gbk = gbk
        self.assemble_unmap = assemble_unmap
        
    def run(self,):
        '''
        description
        '''
        fq = FASTQ_Stats(self.FASTQ_R1, self.FASTQ_R2)
        fq.run()
        sample_name = self.sample_name
        reference = self.reference
        samfile = f'{sample_name}.sam'
        fixmate_bamfile = f'{sample_name}_fixmate.bam'
        pos_srt_bamfile = f'{sample_name}_pos_srt.bam'
        nodup_bamfile = f'{sample_name}_nodup.bam'
        unfiltered_hapall = f'{sample_name}_unfiltered_hapall.vcf'
        mapfix_hapall = f'{sample_name}_mapfix_hapall.vcf'
        filtered_hapall = f'{sample_name}_filtered_hapall.vcf'
        unmapped_read1 = f'{sample_name}_unmapped_R1.fastq'
        unmapped_read2 = f'{sample_name}_unmapped_R2.fastq'
        unmapped_read = f'{sample_name}_unmapped.fastq'
        alignment_vcf_run_summary = []
        os.system(f'samtools faidx {reference}')
        # os.system(f'picard CreateSequenceDictionary REFERENCE={reference} OUTPUT={reference.rsplit(".", 1)[0]}.dict 2> /dev/null')
        os.system(f'bwa index {reference} 2> /dev/null') # Needed for freebayes-parallel.
        if self.nanopore: 
            self.aligner = "Minimap2"
            run_set = f'minimap2 -a -x map-ont -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {reference} {self.FASTQ_R1} -o {samfile}'
            os.system(run_set)
        elif self.paired:
            self.aligner = "BWA"
            # run_set = f'minimap2 -t 8 -a -x sr {reference} {self.FASTQ_R1} {self.FASTQ_R2} -o {samfile}'
            run_set = f'bwa mem -M -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {reference} {self.FASTQ_R1} {self.FASTQ_R2} > {samfile}'
            os.system(run_set)
        else:
            self.aligner = "BWA"
            # run_set = f'minimap2 -t 8 -a -x sr {reference} {self.FASTQ_R1} -o {samfile}'
            run_set = f'bwa mem -M -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {reference} {self.FASTQ_R1} > {samfile}'
            os.system(run_set)
        alignment_vcf_run_summary.append(f'SYSTEM CALL: {run_set} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')

        # http://www.htslib.org/workflow/fastq.html
        samfixmate = f'samtools fixmate -O bam,level=1 -m {samfile} {fixmate_bamfile}'
        samsort = f'samtools sort -l 1 -@8 -o {pos_srt_bamfile} {fixmate_bamfile}'
        sammarkup = f'samtools markdup -f markduplicate_stats.txt -r -O bam,level=1 {pos_srt_bamfile} {nodup_bamfile}'
        os.system(samfixmate)
        os.system(samsort)
        os.system(sammarkup)
        alignment_vcf_run_summary.append(f'SYSTEM CALL: {samfixmate} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
        alignment_vcf_run_summary.append(f'SYSTEM CALL: {samsort} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
        alignment_vcf_run_summary.append(f'SYSTEM CALL: {sammarkup} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')

        # os.remove(fixmate_bamfile)
        # os.remove(pos_srt_bamfile)
        # os.remove(nodup_bamfile)
        # os.system(f'samtools view -Sb {samfile} -o {fixmate_bamfile}')
        # os.system(f'samtools sort {fixmate_bamfile} -o {pos_srt_bamfile}')
        # os.system(f'samtools index {pos_srt_bamfile}')
        # os.system(f'picard MarkDuplicates INPUT={pos_srt_bamfile} OUTPUT={nodup_bamfile} ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv 2> /dev/null')
        # os.system(f'samtools index {nodup_bamfile}')

        # http://www.htslib.org/doc/samtools-markdup.html
        alignment_vcf_run_summary.append(f'NOTE: Read stats gathered by markduplicate_stats.txt -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
        mdf = pd.read_csv('markduplicate_stats.txt', delimiter=':', index_col=0, header=None)
        mdict = mdf.to_dict(orient='dict')[1]
        self.MARKDUPLICATE_COMMAND = mdict['COMMAND']
        self.READS_READ = int(mdict['READ'])
        self.READS_WRITTEN = int(mdict['WRITTEN'])
        self.READS_EXCLUDED = int(mdict['EXCLUDED'])
        self.READS_EXAMINED = int(mdict['EXAMINED'])
        self.READS_PAIRED = int(mdict['PAIRED'])
        self.READS_SINGLE = int(mdict['SINGLE'])
        self.DUPLICATE_PAIR = int(mdict['DUPLICATE PAIR'])
        self.DUPLICATE_SINGLE = int(mdict['DUPLICATE SINGLE'])
        self.DUPLICATE_PAIR_OPTICAL = int(mdict['DUPLICATE PAIR OPTICAL'])
        self.DUPLICATE_SINGLE_OPTICAL = int(mdict['DUPLICATE SINGLE OPTICAL'])
        self.DUPLICATE_NON_PRIMARY = int(mdict['DUPLICATE NON PRIMARY'])
        self.DUPLICATE_NON_PRIMARY_OPTICAL = int(mdict['DUPLICATE NON PRIMARY OPTICAL'])
        self.DUPLICATE_PRIMARY_TOTAL = int(mdict['DUPLICATE PRIMARY TOTAL'])
        self.DUPLICATE_TOTAL = int(mdict['DUPLICATE TOTAL'])
        self.ESTIMATED_LIBRARY_SIZE = int(mdict['ESTIMATED_LIBRARY_SIZE'])
        try:
            self.DUPLICATION_RATIO = float(self.DUPLICATE_TOTAL/self.READS_EXAMINED)
        except ZeroDivisionError:
            self.DUPLICATION_RATIO = 0

        os.system(f'samtools index {nodup_bamfile}')
        if self.nanopore:
            def qual_value_update(vcf):
                df = pd.read_csv(vcf, sep='\t', header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"], comment='#')
                alignment_vcf_run_summary.append(f'NOTE: Nanopore QUAL values increased by 100 to obtain closer values seen with Illumina reads, and allowing VCF files from both platforms to be ran together. -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
                df['QUAL'] = df['QUAL'].apply(lambda x: x + 100)
                header_out = open('v_header.csv', 'w+')
                with open(vcf) as fff:
                    for line in fff:
                        if re.search('^#', line):
                            print(line.strip(), file=header_out)
                header_out.close()
                df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)
                cat_files = ['v_header.csv', 'v_annotated_body.csv']
                name = vcf.replace('.vcf', '')
                qp100 = f'{name}_qp100.vcf'
                with open(qp100, "wb") as outfile:
                    for cf in cat_files:
                        with open(cf, "rb") as infile:
                            outfile.write(infile.read())
                os.remove('v_header.csv')
                os.remove('v_annotated_body.csv')
                return(qp100)

            filtered_hapall = f'{sample_name}_filtered_hapall_nanopore.vcf'
            bcftools_mpileup = f'bcftools mpileup --threads 16 -Ou -f {reference} {nodup_bamfile} | bcftools call --threads 16 -mv -v -Ov -o {unfiltered_hapall}'
            vcffilter = f'vcffilter -f "QUAL > 20" {unfiltered_hapall} > temp1.vcf'
            os.system(bcftools_mpileup)
            os.system(vcffilter)
            alignment_vcf_run_summary.append(f'NOTE: Nanopore - bcftools mpileup used to call SNPs and make VCF files *** -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {bcftools_mpileup} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {vcffilter} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            os.system(f'vcftools --vcf temp1.vcf --remove-indels --recode --recode-INFO-all --out temp2')
            qp100 = qual_value_update('temp2.recode.vcf')
            os.rename(qp100, filtered_hapall)
            os.remove('temp1.vcf')
            os.remove('temp2.recode.vcf')
            # os.remove('temp2.log')
        else:
            chrom_ranges = open("chrom_ranges.txt", 'w')
            for record in SeqIO.parse(reference, "fasta"):
                chrom = record.id
                total_len = len(record.seq)
                min_number = 0
                step = 100000
                if step < total_len:
                    for chunk in range(min_number, total_len, step)[1:]:
                        print("{}:{}-{}".format(chrom, min_number, chunk), file=chrom_ranges)
                        min_number = chunk
                print("{}:{}-{}".format(chrom, min_number, total_len), file=chrom_ranges)
            chrom_ranges.close()
            freebayes_parallel = f'freebayes-parallel chrom_ranges.txt 8 -E -1 -e 1 -u --strict-vcf -f {reference} {nodup_bamfile} > {unfiltered_hapall}' #2>/dev/null'
            os.system(freebayes_parallel)
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {freebayes_parallel} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            # os.system(f'bcftools mpileup --threads 8 -Ou -f {reference} {nodup_bamfile} | bcftools call --threads 8 -mv -Ov -o {unfiltered_hapall}')
            write_fix = open(mapfix_hapall, 'w+')
            with open(unfiltered_hapall, 'r') as unfiltered:
                for line in unfiltered:
                    line = line.strip()
                    new_line = re.sub(r';MQM=', r';MQ=', line)
                    new_line = re.sub(r'ID=MQM,', r'ID=MQ,', new_line)
                    print(new_line, file=write_fix)
                write_fix.close()
            # remove clearly poor positions
            vcffilter_20 = f'vcffilter -f "QUAL > 20" {mapfix_hapall} > {filtered_hapall}'
            os.system(f'vcffilter -f "QUAL > 20" {mapfix_hapall} > {filtered_hapall}')
            alignment_vcf_run_summary.append(f'NOTE: Freebayes MQM= changed to MQ= to stay in the same naming as GATK -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {vcffilter_20} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')

        unmapped_dir = 'unmapped_reads'
        if not os.path.exists(unmapped_dir):
            os.makedirs(unmapped_dir)
        if self.paired and not self.nanopore:
            samtools_fastq_unmapped = f'samtools fastq -f4 -1 {unmapped_read1} -2 {unmapped_read2} --reference {reference} --threads 8 {nodup_bamfile} 2> /dev/null'
            os.system(samtools_fastq_unmapped)
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {samtools_fastq_unmapped} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            os.system(f'pigz {unmapped_read1}')
            os.system(f'pigz {unmapped_read2}')
            shutil.move(f'{unmapped_read1}.gz', unmapped_dir)
            shutil.move(f'{unmapped_read2}.gz', unmapped_dir)
            self.unmapped_read1 = f'{self.cwd}/{unmapped_dir}/{unmapped_read1}.gz'
            self.unmapped_read2 = f'{self.cwd}/{unmapped_dir}/{unmapped_read2}.gz'
            self.unmapped_read_list = [self.unmapped_read1, self.unmapped_read2]
            if self.assemble_unmap:
                assemble = Assemble(FASTQ_R1=self.unmapped_read1, FASTQ_R2=self.unmapped_read2, debug=True)
        elif not self.nanopore:
            samtools_fastq_unmapped = f'samtools fastq -f4 -0 {unmapped_read} --reference {reference} --threads 8 {nodup_bamfile} 2> /dev/null'
            os.system(samtools_fastq_unmapped)
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {samtools_fastq_unmapped} -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            os.system(f'pigz {unmapped_read}')
            shutil.move(f'{unmapped_read}.gz', unmapped_dir)
            self.unmapped_read = f'{self.cwd}/{unmapped_dir}/{unmapped_read}.gz'
            self.unmapped_read_list = [self.unmapped_read]
            if self.assemble_unmap:
                assemble = Assemble(FASTQ_R1=self.unmapped_read, debug=True)

        if not self.assemble_unmap or self.nanopore:
            assembly_message = "skipped assembly"
            alignment_vcf_run_summary.append(f'NOTE: Skipped unmapped read assembly -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            assembly_error = f'{self.cwd}/{unmapped_dir}/skipped_assembly'
            self.unmapped_assemble = None
        else:
            assemble.run()
            try:
                shutil.move(assemble.FASTA, f'{assemble.fastq_name}_unmapped.fasta')
                shutil.move(f'{assemble.fastq_name}_unmapped.fasta', unmapped_dir)
                assemble.FASTA = f'{self.cwd}/{unmapped_dir}/{assemble.fastq_name}_unmapped.fasta'
                assemble.stats(assemble.FASTA)
                self.assemble = assemble
                self.unmapped_assemble = assemble.FASTA
                if not self.debug:
                    shutil.rmtree('spades_assembly')
                assembly_message = f'{assemble.contig_count:,}'
                alignment_vcf_run_summary.append(f'NOTE: SPAdes unmapped read assembly completed -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
            except TypeError:
                self.unmapped_assemble = None
                assembly_error = f'{self.cwd}/{unmapped_dir}/failed_assembly'
                with open(assembly_error, "w") as opened_file:
                    print("see unmapped FASTQ files for troubleshooting", file=opened_file)
                try:
                    shutil.move('spades_assembly/spades.log', unmapped_dir)
                except FileNotFoundError:
                    pass
                if not self.debug and self.assemble_unmap:
                    shutil.rmtree('spades_assembly')
                assembly_message = 'failed assembly'
                alignment_vcf_run_summary.append(f'NOTE: SPAdes unmapped read assembly failed -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
        self.assembly_message = assembly_message

        if self.gbk:
            self.vcf_annotation = VCF_Annotation(gbk_list=self.gbk, vcf_file=filtered_hapall)
            alignment_vcf_run_summary.append(f'IMPORT: VCF_Annotation(gbk_list=self.gbk, vcf_file=filtered_hapall) -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')

        self.zero_coverage = Zero_Coverage(FASTA=reference, bam=nodup_bamfile, vcf=filtered_hapall,)
        alignment_vcf_run_summary.append(f'IMPORT: Zero_Coverage(FASTA=reference, bam=nodup_bamfile, vcf=filtered_hapall,) -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')
        if self.nanopore:
            nanopore_zc_name = re.sub('_zc.vcf', '_nanopore_zc.vcf', self.zero_coverage.zero_coverage_vcf)
            os.rename(self.zero_coverage.zero_coverage_vcf, nanopore_zc_name)
            self.zero_coverage.zero_coverage_vcf = nanopore_zc_name

        alignment = f'alignment_{self.reference_name}'
        if not os.path.exists(alignment):
            os.makedirs(alignment)
        if os.path.exists(unmapped_dir):
            shutil.move(unmapped_dir, alignment)
        files_grab = []
        for files in ('*_nodup.bam', '*_zc.vcf', '*_nodup.bam.bai', '*_annotated.vcf'):
            files_grab.extend(glob.glob(files))
        for each in files_grab:
            shutil.move(each, alignment)
        self.zero_coverage_vcf_file_path = f'{self.cwd}/{alignment}/{os.path.basename(self.zero_coverage.zero_coverage_vcf)}'
        shutil.move(reference, alignment)
        try:
            shutil.move(f'{reference}.fai', alignment)
        except FileNotFoundError:
            pass
        self.reference = f'{self.cwd}/{alignment}/{reference}'
        self.nodup_bamfile = f'{self.cwd}/{alignment}/{nodup_bamfile}'
        self.zc_vcf = self.zero_coverage.zero_coverage_vcf
        if self.gbk:
            self.annotated_vcf = f'{self.cwd}/{alignment}/{os.path.basename(self.vcf_annotation.vcf)}'

        temp_dir = './temp'
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        files_grab = []
        for files in ('*_unmapped*.fastq.gz', '*_all.bam', '*_fixmate.bam', '*_pos_srt.bam', 'markduplicate_stats.txt', '*.bai', '*_filtered_hapall.vcf', '*_mapfix_hapall.vcf', '*_unfiltered_hapall.vcf', '*_filtered_hapall_nanopore.vcf', '*.sam', '*.amb', '*.ann', '*.bwt', '*.pac', '*.fasta.sa', '*_sorted.bam', '*.dict', 'chrom_ranges.txt', '*.fai', 'dup_metrics.csv'):
            files_grab.extend(glob.glob(files))
        for each in files_grab:
            shutil.move(each, temp_dir)

        if self.debug is False:
            shutil.rmtree(temp_dir)
            alignment_vcf_run_summary.append(f'NOTE: Files moved to temp_dir and removed: *_unmapped*.fastq.gz, *_all.bam, *_fixmate.bam, *_pos_srt.bam, markduplicate_stats.txt, *.bai, *_filtered_hapall.vcf, *_mapfix_hapall.vcf, *_unfiltered_hapall.vcf, *_filtered_hapall_nanopore.vcf, *.sam, *.amb, *.ann, *.bwt, *.pac, *.fasta.sa, *_sorted.bam, *.dict, chrom_ranges.txt, *.fai, dup_metrics.csv -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')

        #get program versions
        programs = []
        try:
            minimap2_version = subprocess.run(['minimap2', "--version"], capture_output=True, text=True).stdout.rstrip()
            programs.append(f'Minimap2: {minimap2_version}')
            for line in subprocess.run(['freebayes'], capture_output=True, text=True).stdout.splitlines():
                if 'version' in line:
                    programs.append(f'Freebayes: {line.rstrip().replace("version:  ", "")}')
            try:
                programs.append(subprocess.run(['samtools', 'version'], capture_output=True, text=True).stdout.splitlines()[0])
                programs.append(subprocess.run(['samtools', 'version'], capture_output=True, text=True).stdout.splitlines()[1])
            except IndexError:
                for line in subprocess.run(['samtools'], capture_output=True, text=True).stderr.splitlines():
                    if 'Version' in line:
                        programs.append(f'Samtools: {line}')
            programs.append(assemble.spades_version)
        except:
            pass
        
        self.programs = programs
        self.alignment_vcf_run_summary = alignment_vcf_run_summary


    def latex(self, tex, groups=None):
        blast_banner = Banner(f'Read Mapping against {self.reference_name} using {self.aligner}')
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        
        print(r'\begin{tabular}{ l | l | l | l | l | l }', file=tex)
        print(r'Mapped Paired Reads & Mapped Single Reads & Unmapped Reads & Unmapped Percent & \multicolumn{2}{l}{Unmapped Assembled Contigs} \\', file=tex)
        print(r'\hline', file=tex) 
        mapped_reads = self.READS_PAIRED + self.READS_SINGLE
        total_reads = mapped_reads + self.READS_EXCLUDED
        self.freq_unmapped_reads = self.READS_EXCLUDED / total_reads
        print(f'{self.READS_PAIRED:,} & {self.READS_SINGLE:,} & {self.READS_EXCLUDED:,} & {(self.freq_unmapped_reads*100):,.1f}' + r'\%' + f' & ' + r'\multicolumn{2}{l}{' + f'{self.assembly_message}' + r' } \\', file=tex)
        print(r'\hline', file=tex)
        print(r'\hline', file=tex)
        
        print(r'Duplicate Paired Reads & Duplicate Single Reads & \multicolumn{4}{l}{Duplicate Percent of Mapped Reads} \\', file=tex)
        print(r'\hline', file=tex)
        print(f'{self.DUPLICATE_PAIR:,} & {self.DUPLICATE_SINGLE:,} & ' + r'\multicolumn{4}{l}{' + f'{(self.DUPLICATION_RATIO*100):,.1f}' + r'\%} \\', file=tex)
        print(r'\hline', file=tex)
        print(r'\hline', file=tex)

        print(f'BAM File & Reference Length & Genome with Coverage & Average Depth & No Coverage Bases & Quality SNPs \\\\', file=tex)
        print(r'\hline', file=tex)
        bam = self.zero_coverage.bam.replace('_', '\_')
        print(f'{bam} & {self.zero_coverage.reference_length:,} & {(self.zero_coverage.genome_coverage*100):,.2f}\% & {self.zero_coverage.ave_coverage:,.1f}X & {self.zero_coverage.total_zero_coverage:,} & {self.zero_coverage.good_snp_count:,} \\\\', file=tex)
        print(r'\hline', file=tex)

        if groups:
            print(r'\hline', file=tex)
            print(r'Group Assignment & \multicolumn{4}{l}{ ' + f'{groups}' + r'} \\', file=tex)
            print(r'\hline', file=tex)

        print(r'\end{tabular}', file=tex)

        print(r'\end{adjustbox}', file=tex)
        print(r'\\', file=tex)
        print(r'\end{table}', file=tex)
    
    def excel(self, excel_dict):
        self.aligner
        excel_dict['Aligner'] = f'{self.aligner}'
        excel_dict['Mapped Paired Reads'] = f'{self.READS_PAIRED:,}'
        excel_dict['Mapped Single Reads'] = f'{self.READS_SINGLE:,}'
        excel_dict['Unmapped Reads'] = f'{self.READS_EXCLUDED:,}'
        excel_dict['Unmapped Percent'] = f'{(self.freq_unmapped_reads*100):,.1f}%'
        excel_dict['Unmapped Assembled Contigs'] = f'{self.assembly_message}'
        excel_dict['Duplicate Paired Reads'] = f'{self.DUPLICATE_PAIR:,}'
        excel_dict['Duplicate Single Reads'] = f'{self.DUPLICATE_SINGLE:,}'
        excel_dict['Duplicate Percent of Mapped Reads'] = f'{(self.DUPLICATION_RATIO*100):,.1f}%'
        excel_dict['BAM/Reference File'] = f'{self.zero_coverage.bam} made with {self.reference_name}'
        excel_dict['Reference Length'] = f'{self.zero_coverage.reference_length:,}'
        excel_dict['Genome with Coverage'] = f'{(self.zero_coverage.genome_coverage*100):,.2f}%'
        excel_dict['Average Depth'] = f'{self.zero_coverage.ave_coverage:,.1f}X'
        excel_dict['No Coverage Bases'] = f'{self.zero_coverage.total_zero_coverage:,}'
        excel_dict['Percent Ref with Zero Coverage'] = f'{self.zero_coverage.percent_ref_with_zero_coverage:,.6f}%'
        excel_dict['Quality SNPs'] = f'{self.zero_coverage.good_snp_count:,}'


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Usage:
    alignment_vcf.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta
    alignment_vcf.py -r1 *fastq.gz -r *fasta
    alignment_vcf.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta -g *gbk

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-n', '--SAMPLE_NAME', action='store', dest='SAMPLE_NAME', required=False, help='Force output files to this sample name')
    parser.add_argument('-r1', '--read1', action='store', dest='FASTQ_R1', required=True, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: R2 Illumina read')
    parser.add_argument('-r', '--reference', nargs='*', dest='FASTA', required=False, help='FASTA file to be used as reference.  Multiple can be specified with wildcard')
    parser.add_argument('-o', '--nanopore', action='store_true', dest='nanopore', default=False, help='Beta, if true run alignment with minimap2 map-ont option')
    parser.add_argument('-b', '--gbk', nargs='*', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file.  Multiple can be specified with wildcard')
    parser.add_argument('-assemble_unmap', '--assemble_unmap', action='store_true', dest='assemble_unmap', help='skip assembly of unmapped reads')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    def concat_fasta(FASTAs):
        basenames=[]
        fasta_name=[]
        for each in FASTAs:
            basenames.append(os.path.basename(each))
            fasta_name.append(re.sub('\..*', '', os.path.basename(each)))
        fastas_used = ", ".join(basenames)
        #concatenate FASTA list
        concatenated_TEMP = f'{"_".join(fasta_name)}.temp'
        concatenated_FASTA = f'{"_".join(fasta_name)}.fasta'
        with open(concatenated_TEMP,'wb') as wfd:
            for each in FASTAs:
                with open(each,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
        shutil.move(concatenated_TEMP, concatenated_FASTA)
        return concatenated_FASTA, fastas_used
    concatenated_FASTA, fastas_used = concat_fasta(args.FASTA) # -f option for FASTA will take wildcard for multiple FASTAs

    alignment = Alignment(SAMPLE_NAME=args.SAMPLE_NAME, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference=concatenated_FASTA, nanopore=args.nanopore, gbk=args.gbk, assemble_unmap=args.assemble_unmap, debug=args.debug)
    alignment.run()

    #Latex report
    latex_report = Latex_Report(alignment.sample_name)
    alignment.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(alignment.sample_name)
    excel_stats.excel_dict["FASTA/s"] = fastas_used
    alignment.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*.aux', '*.log', '*tex', '*png', '*out', '*_all.bam', '*.bai', '*_filtered_hapall.vcf', '*_mapfix_hapall.vcf', '*_unfiltered_hapall.vcf', '*.sam', '*.amb', '*.ann', '*.bwt', '*.pac', '*.fasta.sa', '*_sorted.bam', '*.dict', 'chrom_ranges.txt', 'dup_metrics.csv', '*.fai'):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        shutil.move(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# Created March 2021 by Tod Stuber
