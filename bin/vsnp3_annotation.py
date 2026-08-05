#!/usr/bin/env python3

from vsnp3_version import __version__

import os
import shutil
import sys
import argparse
import textwrap
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Data import CodonTable, IUPACData
from collections import ChainMap
from Bio.SeqFeature import FeatureLocation, CompoundLocation

# One letter to three letter amino acid names, for HGVS protein notation.
_AA_1TO3 = dict(IUPACData.protein_letters_1to3)
_AA_1TO3['_'] = 'Ter'
_AA_1TO3['*'] = 'Ter'


class Annotation():
    '''
    Locate a genomic position within an annotated feature and report the amino acid
    consequence of a single nucleotide substitution.

    Codons are taken from the feature's own coding sequence via feature.extract(),
    which Biopython returns spliced and in gene orientation.  Nothing here re-derives
    strand, exon boundaries or reading frame by hand, because every historical defect
    in this module came from doing exactly that.

    Conventions kept for backward compatibility with the SNP tables:
      - stop codons are written '_' , not Biopython's '*'
      - a position outside every annotated feature leaves gene = 'unlisted gene'
    '''

    # Feature types this module will locate a position within.
    ANNOTATED_TYPES = ("CDS", "tRNA", "rRNA", "repeat_region", "mobile_element", "ncRNA")

    # Reference validation is per gbk file, not per Annotation instance.  step2 builds
    # one Annotation per VCF, and revalidating a 4000 CDS reference each time is waste.
    _validated = {}

    def __init__(self, gbk_list=None, validate=True):

        self.aa_code = { \
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', \
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', \
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', \
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', \
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', \
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', \
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', \
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', \
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', \
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', \
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', \
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', \
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', \
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', \
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', \
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

        cwd = os.getcwd()
        gbk_dict_list=[]
        self._feature_cache = {}   # (chrom, id(feature)) -> per feature coding info
        # Chroms seen in the data that no supplied gbk describes.  Reported once at
        # the end of a run rather than terminating it at the first occurrence.
        self.chroms_without_records = set()

        # Ensure gbk_list is not None and handle file processing safely
        if gbk_list is None:
            gbk_list = []

        # Sorted so that when two gbk files define the same chrom id the winner is
        # the same on every machine.  ChainMap gives first-wins, and the list used to
        # arrive in directory order.
        source_of_chrom = {}
        for each in sorted(gbk_list):
            try: #IF gbk provided as path variable copy local
                records = SeqIO.to_dict(SeqIO.parse(each, "genbank"))
            except (FileNotFoundError, PermissionError, OSError) as e:
                print(f"Warning: Could not process file {each}: {e}")
                continue
            except shutil.SameFileError:
                continue
            for chrom in records:
                if chrom in source_of_chrom:
                    # Sorting makes this deterministic but not correct: the reference
                    # directory is ambiguous and only its owner can resolve that.
                    print(f"Warning: chrom {chrom} is defined in both "
                          f"{os.path.basename(source_of_chrom[chrom])} and "
                          f"{os.path.basename(each)}; annotation will use "
                          f"{os.path.basename(source_of_chrom[chrom])}")
                else:
                    source_of_chrom[chrom] = each
            gbk_dict_list.append(records)
        self.gbk_dict = dict(ChainMap(*gbk_dict_list)) #merge dictionaries

        if validate:
            for each in gbk_list:
                self.validate_reference(each)
            self.self_test()

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _is_coding(feature):
        '''A codon is only meaningful for a CDS that is not a pseudogene.'''
        return (feature.type == "CDS"
                and 'pseudo' not in feature.qualifiers
                and 'pseudogene' not in feature.qualifiers)

    @staticmethod
    def _start_codons(table):
        try:
            return CodonTable.unambiguous_dna_by_id[int(table)].start_codons
        except (KeyError, ValueError):
            return CodonTable.unambiguous_dna_by_id[11].start_codons

    def _translate_codon(self, codon, table, is_first_codon=False):
        '''
        Translate one codon, honouring the table's initiator rule.  An alternative
        start codon such as GTG or TTG reads as M at residue 1; 38% of H37Rv CDS
        start with something other than ATG, so this is not an edge case.
        '''
        codon = str(codon).upper()
        if len(codon) != 3:
            return 'unfound_ref_AA'
        if is_first_codon and codon in self._start_codons(table):
            return 'M'
        try:
            aa = str(Seq(codon).translate(table=table))
        except Exception:
            return 'ambiguous'
        if aa == '*':
            return '_'          # keep the existing SNP table convention
        if aa == 'X':
            return 'ambiguous'
        return aa

    def _feature_info(self, chrom, feature):
        '''
        Per feature coding information, built once and cached.

        offsets  genomic 0-based position -> offset from the feature's 5' end.
                 list(feature.location) is already in gene order and splice aware,
                 and is reversed for minus strand features, so this needs no
                 strand or exon arithmetic of its own.
        '''
        key = (chrom, id(feature))
        info = self._feature_cache.get(key)
        if info is None:
            record = self.gbk_dict[chrom]
            offsets = {int(pos): i for i, pos in enumerate(feature.location)}
            cds_seq = str(feature.extract(record.seq))
            table = feature.qualifiers.get('transl_table', ['11'])[0]
            codon_start = int(feature.qualifiers.get('codon_start', ['1'])[0]) - 1
            protein = feature.qualifiers.get('translation', [None])[0]
            info = (offsets, cds_seq, protein, table, codon_start)
            self._feature_cache[key] = info
        return info

    # ------------------------------------------------------- reference validation

    def validate_reference(self, gbk_path):
        '''
        Translate every CDS carrying a /translation and compare with the qualifier.

        On a sound reference this agrees for every gene, so any disagreement means the
        reference itself is damaged or uses a translation table Biopython was not told
        about.  This checks the reference only - see self_test() for the check that
        would have caught the 3.34/3.35 codon regression in this module.
        '''
        try:
            stamp = (os.path.abspath(gbk_path), os.path.getmtime(gbk_path))
        except OSError:
            return
        if stamp in Annotation._validated:
            return

        checked = mismatched = uncheckable = no_translation = no_product = total = 0
        try:
            for record in SeqIO.parse(gbk_path, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    total += 1
                    if 'product' not in feature.qualifiers:
                        no_product += 1
                    translation = feature.qualifiers.get('translation', [None])[0]
                    if translation is None:
                        no_translation += 1
                        continue
                    table = feature.qualifiers.get('transl_table', ['11'])[0]
                    try:
                        computed = str(feature.extract(record.seq).translate(table=table, cds=True))
                    except Exception:
                        uncheckable += 1
                        continue
                    checked += 1
                    if computed != translation:
                        mismatched += 1
        except (FileNotFoundError, PermissionError, OSError, ValueError):
            return

        Annotation._validated[stamp] = True
        name = os.path.basename(gbk_path)
        if total == 0:
            print(f"Warning: {name} contains no CDS features; no amino acid annotation is possible")
            return
        if mismatched:
            print(f"### ANNOTATION SELF CHECK FAILED: {name}")
            print(f"###   {mismatched} of {checked} CDS do not match their own /translation qualifier.")
            print(f"###   Amino acid calls from this reference should not be trusted.")
        if no_translation:
            print(f"Note: {name} has {no_translation} of {total} CDS with no /translation; "
                  f"reference amino acids for those are computed, not verified")
        if no_product:
            print(f"Note: {name} has {no_product} of {total} CDS with no /product; "
                  f"the product column falls back to /locus_tag")

    def _codon_call(self, chrom, feature, position_0based, snp_nt):
        """
        Coding consequence of a substitution for ONE feature.  Returns a dict, or
        None when the position has no codon in this feature.

        Sets nothing on self, so it can be reused for every overlapping CDS.
        """
        offsets, cds_seq, protein, table, codon_start = self._feature_info(chrom, feature)
        gene_offset = offsets.get(position_0based)
        if gene_offset is None:
            return None
        gene_offset -= codon_start
        if gene_offset < 0:                      # ahead of the coding start
            return None
        codon_number, nt_index = divmod(gene_offset, 3)
        left = codon_start + codon_number * 3
        codon = cds_seq[left:left + 3]
        if len(codon) < 3:                       # ragged 3' end of a partial CDS
            return None

        # The table's initiator rule only applies if the reference codon really is an
        # initiator.  Otherwise a variant that happened to create GTG or TTG would be
        # reported as M at a position that never started a protein.
        ref_is_initiator = (codon_number == 0
                            and codon.upper() in self._start_codons(table))

        # /translation is the authoritative reference protein when the record carries
        # one.  Translating the codon is only a fallback for references built without.
        if protein is not None and codon_number < len(protein):
            ref_aa = '_' if protein[codon_number] == '*' else protein[codon_number]
        else:
            ref_aa = self._translate_codon(codon, table, ref_is_initiator)

        call = dict(gene=(feature.qualifiers.get('gene', [None])[0]
                          or feature.qualifiers.get('locus_tag', [None])[0]
                          or 'unlisted gene'),
                    codon=codon, codon_number=codon_number, nt_index=nt_index,
                    gene_nt_pos=gene_offset + 1, ref_aa=ref_aa,
                    start_codon=(codon_number == 0),
                    snp_base_code="SNP nt not provided", snp_aa="SNP nt not provided",
                    mutation_type="n/a", hgvs_c="n/a", hgvs_p="n/a")

        if snp_nt is None or snp_nt in ('', '.'):
            return call
        if len(str(snp_nt)) != 1:
            # Insertions, deletions and complex alleles are outside what a single
            # codon substitution can express.  Say so instead of returning
            # 'ambiguous' as a side effect of a failed dictionary lookup.
            call.update(snp_base_code=str(snp_nt), snp_aa="n/a",
                        mutation_type="indel, not evaluated")
            return call

        # On a minus strand gene the VCF ALT is given on the plus strand and has to be
        # complemented before it goes into a gene oriented codon.
        substitution = str(snp_nt).upper()
        if feature.location.strand == -1:
            substitution = str(Seq(substitution).complement())
        codon_list = list(codon)
        codon_list[nt_index] = substitution
        snp_base_code = "".join(codon_list)
        snp_aa = self._translate_codon(snp_base_code, table, ref_is_initiator)

        # Residue 1 is a special case: a variant that swaps one valid initiator for
        # another is genuinely silent, anything else destroys the start.
        if snp_aa == "ambiguous":
            mutation_type = "unsure-ambiguous"
        elif ref_aa == "unfound_ref_AA":
            mutation_type = ""
        elif ref_is_initiator and snp_aa != ref_aa:
            mutation_type = "start codon lost"
        elif ref_aa == snp_aa:
            mutation_type = "silent mutation"
        else:
            mutation_type = "nonsynonymous"

        call.update(snp_base_code=snp_base_code, snp_aa=snp_aa,
                    mutation_type=mutation_type)

        # HGVS, the form the WHO catalogue and TB-Profiler use.  Deliberately free of
        # whitespace and '=' because this is written into the VCF ID field: the strict
        # synonymous shorthand 'p.Ser315=' would break every key=value parser.
        call['hgvs_c'] = f"c.{gene_offset + 1}{codon[nt_index]}>{substitution}"
        ref3 = _AA_1TO3.get(ref_aa)
        snp3 = _AA_1TO3.get(snp_aa)
        if ref3 and snp3:
            call['hgvs_p'] = f"p.{ref3}{codon_number + 1}{snp3}"
        return call

    def self_test(self, sample=60):
        '''
        Exercise run() against the reference's own /translation and confirm that the
        codon this module reports actually translates to the residue the reference
        says is there.

        validate_reference above checks the reference; this checks the code.  The two
        are different failures and only this one would have caught the minus-strand
        codon defect: in 3.34/3.35 the reported codon was the plus-strand reading, so
        translating it disagrees with /translation on every minus-strand gene.
        '''
        failures = []
        tested = 0
        for chrom, record in self.gbk_dict.items():
            candidates = [f for f in record.features
                          if self._is_coding(f) and 'translation' in f.qualifiers]
            if not candidates:
                continue
            # deterministic spread rather than random, so the check is reproducible
            step = max(1, len(candidates) // max(1, sample))
            for feature in candidates[::step][:sample]:
                protein = feature.qualifiers['translation'][0]
                table = feature.qualifiers.get('transl_table', ['11'])[0]
                order = [int(p) for p in feature.location]
                if len(order) < 6:
                    continue
                # a position in the middle of the gene, away from start/stop
                offset = (len(order) // 6) * 3
                codon_number = offset // 3
                if codon_number >= len(protein):
                    continue
                self.run(f"{chrom}:{order[offset] + 1}", None)
                if not self.is_coding or self.overlapping_features > 1:
                    continue
                if self.aa_residue_pos != codon_number + 1:
                    continue          # a different overlapping feature was chosen
                tested += 1
                expected = protein[codon_number]
                if expected == '*':
                    expected = '_'
                observed = self._translate_codon(self.reference_base_code, table)
                if observed != expected:
                    failures.append((chrom, order[offset] + 1,
                                     feature.qualifiers.get('gene', ['?'])[0],
                                     self.reference_base_code, observed, expected,
                                     self.direction))
        if failures:
            print("### ANNOTATION SELF TEST FAILED "
                  f"({len(failures)} of {tested} sampled codons)")
            for chrom, pos, gene, codon, observed, expected, direction in failures[:5]:
                print(f"###   {chrom}:{pos} {gene} ({direction}): reported codon {codon} "
                      f"translates to {observed}, reference protein has {expected}")
            reverse = sum(1 for f in failures if f[6] == "reverse gene")
            print(f"###   {reverse} of {len(failures)} are on the reverse strand")
            print("###   Amino acid columns in this run are not trustworthy.")
        return not failures

    # -------------------------------------------------------------------- run

    def run(self, abs_pos, snp_nt):
        '''
        When having both abs_pos and snp_nt, the snp_nt allows updating the amino acid to the mutant call
        # abs_pos='NC_017250.1:264518', snp_nt='T'
        # abs_pos='NC_017251.1:396173', snp_nt='T'
        '''

        # Validate input parameters
        if not abs_pos or ':' not in abs_pos:
            print("Error: Invalid abs_pos format. Expected format: 'chrom:position'")
            return

        try:
            chrom, position = abs_pos.split(':', 1)  # Split only on first ':'
            position = int(position)
        except ValueError:
            print(f"Error: Could not parse position from abs_pos: {abs_pos}")
            return

        self.chrom = chrom
        self.position = position
        self.cds_nt_start = "n/a"
        self.cds_nt_end = "n/a"
        self.aa_residue_pos = "n/a"
        self.mutation_type = "n/a"
        self.snp_nt = snp_nt
        self.reference_base_code = "n/a"
        self.snp_base_code = "SNP nt not provided"
        self.ref_aa = "n/a"
        self.snp_aa = "n/a"
        self.gene = "unlisted gene"
        self.product = "No annotated product"
        self.aa_pos = "n/a"
        self.feature_found = False
        self.direction = "n/a"
        self.feature_type = "n/a"
        self.is_coding = False
        self.overlapping_features = 0
        self.overlapping_genes = "n/a"
        self.cds_overlap = "n/a"
        self.start_codon = False
        self.gene_nt_pos = "n/a"
        self.hgvs_c = "n/a"
        self.hgvs_p = "n/a"
        self.all_consequences = "n/a"

        # Convert 1-based genomic position to 0-based for BioPython compatibility
        position_0based = position - 1

        try:
            if chrom not in self.gbk_dict:
                raise KeyError(f"Chromosome {chrom} not found in provided GenBank files")

            # Collect every annotated feature containing the position.  Comparing
            # 0-based against a half open [start, end) interval; the old
            # 'position in range(part.start, part.end)' test compared the 1-based
            # position and so lost the final base of every feature while wrongly
            # claiming the base immediately upstream of it.
            matches = []
            for feature in self.gbk_dict[chrom].features:
                if feature.type not in self.ANNOTATED_TYPES:
                    continue
                location_parts = list(feature.location.parts) or [feature.location]
                for part_index, part in enumerate(location_parts):
                    if int(part.start) <= position_0based < int(part.end):
                        matches.append((feature, part_index, part))
                        break

            if not matches:
                return

            # Prefer a real CDS over an overlapping non-coding or pseudogene feature.
            # The old code returned on whichever feature appeared first in the file,
            # so a repeat_region could mask the gene the position actually sits in.
            # Deterministic primary: longest coding sequence, then lowest start, then
            # locus_tag.  File order is deliberately NOT used - the Mallard and Owl
            # references list PA/PA-X and NS1/NEP in opposite orders, so file order
            # answers the same position two different ways for the same organism.
            # Longest agrees across both and picks the canonical product each time
            # (PB1 over PB1-F2, PA over PA-X, M1 over M2, NS1 over NEP).
            def _primary_rank(match):
                candidate = match[0]
                coding_length = sum(int(part.end) - int(part.start)
                                    for part in candidate.location.parts)
                return (-coding_length, int(candidate.location.start),
                        candidate.qualifiers.get('locus_tag', [''])[0])

            coding_matches = sorted((m for m in matches if self._is_coding(m[0])),
                                    key=_primary_rank)
            feature, part_index, part = (coding_matches or matches)[0]

            # Preferring a CDS over a non-coding feature is defensible.  Choosing
            # between two overlapping CDS is not - there is no single correct answer,
            # so report every candidate rather than silently picking one.
            if len(matches) > 1:
                # Coding features first in primary order, then the rest, so the list is
                # identical for two references that annotate the same position with the
                # features listed in a different order.
                ordered = coding_matches + [m for m in matches if m not in coding_matches]
                names = []
                for other, _, _ in ordered:
                    label = (other.qualifiers.get('gene', [None])[0]
                             or other.qualifiers.get('locus_tag', [None])[0]
                             or other.type)
                    if label not in names:
                        names.append(label)
                self.overlapping_genes = ",".join(names)
            if len(coding_matches) > 1:
                # Kept token-safe: this is written into the VCF ID field.  The
                # gene names themselves are in overlapping_genes.
                self.cds_overlap = "yes"

            self.feature_found = True
            self.overlapping_features = len(matches)
            self.feature_type = feature.type
            self.is_coding = self._is_coding(feature)
            self.direction = "reverse gene" if feature.location.strand == -1 else "forward gene"
            # Report the bounds of the whole gene, not of the matching exon.
            self.cds_nt_start = int(feature.location.start)
            self.cds_nt_end = int(feature.location.end)

            if 'gene' in feature.qualifiers:
                self.gene = feature.qualifiers['gene'][0]

            # Product, with a fallback chain that is actually reachable.  This was
            # previously written as nested try/except guarded by an 'in' test, so an
            # absent qualifier raised nothing and the fallbacks were dead code.
            if 'product' in feature.qualifiers:
                self.product = feature.qualifiers['product'][0]
            elif 'locus_tag' in feature.qualifiers:
                self.product = feature.qualifiers['locus_tag'][0]
            elif 'label' in feature.qualifiers:
                self.product = f'{feature.type}, {feature.qualifiers["label"][0]}'
            else:
                self.product = f'{feature.type}, product_unknown'

            # A codon is only meaningful inside a coding feature.  tRNA, rRNA,
            # ncRNA, repeat_region, mobile_element and pseudogenes previously came
            # back with a confident silent/nonsynonymous call.
            if not self.is_coding:
                self.mutation_type = "non-coding"
                self.snp_base_code = "n/a"
                return

            call = self._codon_call(chrom, feature, position_0based, snp_nt)
            if call is None:
                self.mutation_type = "non-coding"
                self.snp_base_code = "n/a"
                return

            self.aa_residue_pos = call['codon_number'] + 1
            self.aa_pos = call['nt_index'] + 1
            self.start_codon = call['start_codon']
            self.reference_base_code = call['codon']
            self.ref_aa = call['ref_aa']
            self.snp_base_code = call['snp_base_code']
            self.snp_aa = call['snp_aa']
            self.mutation_type = call['mutation_type']
            self.gene_nt_pos = call['gene_nt_pos']
            self.hgvs_c = call['hgvs_c']
            self.hgvs_p = call['hgvs_p']

            # Each overlapping CDS is a separate true statement about this position -
            # in influenza a SNP in the PA/PA-X or M1/M2 overlap really does change two
            # proteins, in different frames.  Choosing a primary is a presentation
            # decision, so report the others rather than discarding them.
            if len(coding_matches) > 1:
                consequences = []
                for other, _, _ in coding_matches:
                    other_call = self._codon_call(chrom, other, position_0based, snp_nt)
                    if other_call is None:
                        continue
                    detail = (other_call['hgvs_p'] if other_call['hgvs_p'] != 'n/a'
                              else other_call['mutation_type'] or 'no call')
                    consequences.append(f"{other_call['gene']}:{detail}")
                if consequences:
                    self.all_consequences = "|".join(consequences)

            return

        except KeyError:
            # A chrom with no GenBank record is a normal condition, not a fatal one:
            # multi-contig references where only some contigs have a gbk, plasmids
            # and phage all reach here.  This used to sys.exit(0), so one such
            # contig terminated the whole of step 2 mid-run with a SUCCESS status -
            # no tables, no trees, and nothing for a scheduler or wrapper to detect.
            # The position is now reported as unannotated and the run continues.
            self.mutation_type = 'no reference record for chrom'
            self.chroms_without_records.add(self.chrom)
            return
        except Exception as e:
            print(f"Unexpected error in annotation processing: {e}")
            return

        return


class VCF_Annotation():
    '''
    '''
    def __init__(self, vcf_file=None,):
        self.vcf_file = vcf_file


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------

    Used by vSNP to annotate VCF files and SNP tables

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-g', '--gbk', nargs='*', dest='gbk', default=None, required=True, help='Full path to .gbk files.  Wildcard can be used')
    parser.add_argument('-p', '--abs_pos', action='store', dest='abs_pos', default=None, required=True, help='Absolute position, example NC_017250.1:264518')
    parser.add_argument('-s', '--snp_nt', action='store', dest='snp_nt', default=None, required=False, help='Position alt call, example T')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    print("\n")

    try:
        annotation = Annotation(gbk_list=args.gbk)
        annotation.run(args.abs_pos, args.snp_nt)

        # Safely print results with getattr fallbacks to maintain exact output format
        abs_pos_safe = getattr(args, 'abs_pos', 'unknown')
        print(abs_pos_safe)

        print(f'cds_nt_start: {getattr(annotation, "cds_nt_start", "n/a")}')
        print(f'cds_nt_end: {getattr(annotation, "cds_nt_end", "n/a")}')
        print(f'gene: {getattr(annotation, "gene", "unlisted gene")}')
        print(f'product: {getattr(annotation, "product", "No annotated product")}')
        print(f'feature_type: {getattr(annotation, "feature_type", "n/a")}')
        print(f'aa_residue_pos: {getattr(annotation, "aa_residue_pos", "n/a")}')
        print(f'snp_nt: {getattr(annotation, "snp_nt", "n/a")}')
        print(f'aa_pos: {getattr(annotation, "aa_pos", "n/a")}')
        print(f'reference base code: {getattr(annotation, "reference_base_code", "n/a")}')
        print(f'snp_base_code: {getattr(annotation, "snp_base_code", "SNP nt not provided")}')
        print(f'ref_aa: {getattr(annotation, "ref_aa", "n/a")}')
        print(f'snp_aa: {getattr(annotation, "snp_aa", "n/a")}')
        print(f'mutation_type: {getattr(annotation, "mutation_type", "n/a")}')
        print(f'Gene direction: {getattr(annotation, "direction", "n/a")}')
        print(f'gene_nt_pos: {getattr(annotation, "gene_nt_pos", "n/a")}')
        print(f'hgvs_c: {getattr(annotation, "hgvs_c", "n/a")}')
        print(f'hgvs_p: {getattr(annotation, "hgvs_p", "n/a")}')
        print(f'all_consequences: {getattr(annotation, "all_consequences", "n/a")}')
        print(f'overlapping_features: {getattr(annotation, "overlapping_features", 0)}')
        print(f'overlapping_genes: {getattr(annotation, "overlapping_genes", "n/a")}')
        if getattr(annotation, "cds_overlap", "n/a") == "yes":
            print(f'cds_overlap: ambiguous, position lies in multiple CDS '
                  f'({getattr(annotation, "overlapping_genes", "")})')
        else:
            print('cds_overlap: n/a')

        # Final mutation notation with safe attribute access
        ref_aa = getattr(annotation, "ref_aa", "n/a")
        aa_residue_pos = getattr(annotation, "aa_residue_pos", "n/a")
        snp_aa = getattr(annotation, "snp_aa", "n/a")
        print(f'{ref_aa}{aa_residue_pos}{snp_aa}')

    except Exception as e:
        print(f"Error: Failed to create or run annotation: {e}")
        sys.exit(1)

# Created 2021 by Tod Stuber
