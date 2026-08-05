#!/usr/bin/env python3
"""
Regression tests for vsnp3_annotation.py.

Sections 1-3 are self-contained and need no data files; they are the ones that would
have caught the minus-strand codon defect that shipped in 3.34/3.35.
Sections 4-6 use the shipped references and the four sample VCFs when present.

    python3 tests/test_annotation.py
    VSNP3_REQUIRE_DATA=1 python3 tests/test_annotation.py   # skips become failures

Sections 4-6 need the reference set and sample VCFs, which are not in this repo.
Point VSNP3_DEPS and VSNP3_VCFS at them.  Without VSNP3_REQUIRE_DATA those
sections skip, and a skip does not fail the run -- so on a machine without the
data this suite reports success while exercising about half of what it claims.
Set VSNP3_REQUIRE_DATA=1 in CI for the tiers whose data is available.
"""

import os
import sys
import glob
import random

sys.path.insert(0, os.environ.get(
    "VSNP3_BIN",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin")))

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from Bio.Data import CodonTable
from vsnp3_annotation import Annotation

START_CODONS = CodonTable.unambiguous_dna_by_id[11].start_codons

# No developer-machine defaults: an absent path must be visibly absent, not a
# path that happens to exist on one laptop.
DEP = os.environ.get("VSNP3_DEPS", "")
VCFDIR = os.environ.get("VSNP3_VCFS", "")
REQUIRE_DATA = os.environ.get("VSNP3_REQUIRE_DATA", "") not in ("", "0", "false")

PASS, FAIL, SKIP = [], [], []


def check(name, condition, detail=""):
    (PASS if condition else FAIL).append(name)
    print(f"  {'PASS' if condition else 'FAIL'}  {name}" + (f"   {detail}" if detail and not condition else ""))
    return bool(condition)


def skip(name, why):
    if REQUIRE_DATA:
        FAIL.append(name)
        print(f"  FAIL  {name}   ({why}; VSNP3_REQUIRE_DATA is set)")
        return
    SKIP.append(name)
    print(f"  SKIP  {name}   ({why})")


def norm(aa):
    return "*" if aa in ("_", "*") else aa


# --------------------------------------------------------------------- section 1

def synthetic_reference(length=600, seed=7):
    random.seed(seed)
    return "".join(random.choice("ATGC") for _ in range(length))


def truth(genome, p0, alt, parts_gene_order, strand):
    """
    Independent ground truth.  parts_gene_order is [(start, end), ...] in 5'->3' GENE
    order, so for a minus strand feature the highest-coordinate exon comes first.
    Build the gene sequence and the gene-ordered position list, then index into them.
    """
    if strand == -1:
        gene = "".join(str(Seq(genome[s:e]).reverse_complement()) for s, e in parts_gene_order)
        order = [p for s, e in parts_gene_order for p in range(e - 1, s - 1, -1)]
        a = str(Seq(alt).complement())
    else:
        gene = "".join(genome[s:e] for s, e in parts_gene_order)
        order = [p for s, e in parts_gene_order for p in range(s, e)]
        a = alt
    off = order.index(p0)
    ci, idx = divmod(off, 3)
    codon = gene[3 * ci:3 * ci + 3]
    if len(codon) < 3:
        return None
    mut = codon[:idx] + a + codon[idx + 1:]
    ref_aa = str(Seq(codon).translate())
    alt_aa = str(Seq(mut).translate())
    # Table 11 initiator rule, applied only when the reference codon is itself an
    # initiator -- the same condition the module uses.
    if ci == 0 and codon in START_CODONS:
        ref_aa = "M"
        if mut in START_CODONS:
            alt_aa = "M"
    return ci + 1, ref_aa, alt_aa


def section1_synthetic():
    print("\n[1] synthetic reference: both strands, all three codon offsets, ragged lengths")
    genome = synthetic_reference()
    cases = [("plus  len%%3==0", 0, 90, 1), ("minus len%%3==0", 150, 240, -1),
             ("minus len%%3==1", 150, 238, -1), ("minus len%%3==2", 150, 239, -1),
             ("plus  len%%3==1", 300, 388, 1)]
    for label, s, e, strand in cases:
        rec = SeqRecord(Seq(genome), id="CHR1", name="CHR1")
        rec.features = [SeqFeature(FeatureLocation(s, e, strand=strand), type="CDS",
                                   qualifiers={"gene": ["g"], "product": ["p"]})]
        ann = Annotation(gbk_list=[], validate=False)
        ann.gbk_dict = {"CHR1": rec}
        n = ok = 0
        first = None
        for p0 in range(s, e):
            for alt in "ATGC":
                if alt == genome[p0]:
                    continue
                t = truth(genome, p0, alt, [(s, e)], strand)
                if t is None:
                    continue
                n += 1
                ann.run(f"CHR1:{p0 + 1}", alt)
                got = (ann.aa_residue_pos, norm(ann.ref_aa), norm(ann.snp_aa))
                syn_ok = ((t[1] == t[2]) == (ann.mutation_type in ("silent mutation",)))
                if got == t and syn_ok:
                    ok += 1
                elif first is None:
                    first = f"pos {p0+1} {genome[p0]}>{alt}: got {got} want {t} type={ann.mutation_type}"
        check(f"{label}  {ok}/{n} exact", ok == n, first or "")


def section2_spliced():
    print("\n[2] spliced CDS (CompoundLocation), the influenza M2/NEP/PA-X shape")
    genome = synthetic_reference(1200, seed=13)
    # exon1 length 26 -> frame shifted continuation, the M2 case
    for label, parts, strand in (("plus  exon1=26 (frame shift)", [(0, 26), (700, 968)], 1),
                                 ("plus  exon1=30 (in frame)",    [(0, 30), (700, 1036)], 1),
                                 ("minus exon1=26 (frame shift)", [(0, 26), (700, 968)], -1)):
        ordered = parts if strand == 1 else list(reversed(parts))
        loc = CompoundLocation([FeatureLocation(s, e, strand=strand) for s, e in ordered])
        rec = SeqRecord(Seq(genome), id="CHR1", name="CHR1")
        rec.features = [SeqFeature(loc, type="CDS", qualifiers={"gene": ["m2"], "product": ["p"]})]
        ann = Annotation(gbk_list=[], validate=False)
        ann.gbk_dict = {"CHR1": rec}
        n = ok = 0
        first = None
        for s, e in ordered:
            for p0 in range(s, e):
                for alt in "ATGC":
                    if alt == genome[p0]:
                        continue
                    t = truth(genome, p0, alt, ordered, strand)
                    if t is None:
                        continue
                    n += 1
                    ann.run(f"CHR1:{p0 + 1}", alt)
                    got = (ann.aa_residue_pos, norm(ann.ref_aa), norm(ann.snp_aa))
                    if got == t:
                        ok += 1
                    elif first is None:
                        first = f"pos {p0+1}: got {got} want {t}"
        check(f"{label}  {ok}/{n} exact", ok == n, first or "")


def section3_boundaries_and_noncoding():
    print("\n[3] feature boundaries, non-coding features, pseudogenes, indels")
    genome = synthetic_reference(600, seed=21)
    rec = SeqRecord(Seq(genome), id="CHR1", name="CHR1")
    rec.features = [
        SeqFeature(FeatureLocation(150, 240, strand=-1), type="CDS",
                   qualifiers={"gene": ["minusGene"], "product": ["mp"]}),
        SeqFeature(FeatureLocation(300, 325, strand=1), type="tRNA",
                   qualifiers={"gene": ["trnA"], "product": ["tRNA-Ala"]}),
        SeqFeature(FeatureLocation(400, 430, strand=1), type="CDS",
                   qualifiers={"gene": ["pseudoG"], "product": ["ps"], "pseudo": [""]}),
        SeqFeature(FeatureLocation(500, 530, strand=1), type="CDS",
                   qualifiers={"gene": ["noProd"], "locus_tag": ["LT_0001"]}),
    ]
    ann = Annotation(gbk_list=[], validate=False)
    ann.gbk_dict = {"CHR1": rec}

    # CDS occupies 0-based [150,240) == 1-based 151..240
    ann.run("CHR1:150", "A")
    check("base upstream of CDS is not claimed", not ann.feature_found,
          f"feature_found={ann.feature_found} residue={ann.aa_residue_pos}")
    ann.run("CHR1:151", "A")
    check("first 1-based base of CDS annotated", ann.feature_found and ann.aa_residue_pos == 30,
          f"residue={ann.aa_residue_pos}")
    ann.run("CHR1:240", "A")
    check("final base of CDS annotated (minus: residue 1)",
          ann.feature_found and ann.aa_residue_pos == 1, f"residue={ann.aa_residue_pos}")
    ann.run("CHR1:241", "A")
    check("base past CDS end not claimed", not ann.feature_found)

    ann.run("CHR1:313", "A")
    check("tRNA reports gene/product but no amino acid",
          ann.feature_found and ann.gene == "trnA" and ann.ref_aa == "n/a"
          and ann.snp_aa == "n/a" and ann.mutation_type == "non-coding",
          f"ref_aa={ann.ref_aa} snp_aa={ann.snp_aa} type={ann.mutation_type}")

    ann.run("CHR1:415", "A")
    check("/pseudo CDS gets no amino acid",
          ann.ref_aa == "n/a" and ann.mutation_type == "non-coding",
          f"ref_aa={ann.ref_aa} type={ann.mutation_type}")

    ann.run("CHR1:515", "A")
    check("product falls back to locus_tag", ann.product == "LT_0001", f"product={ann.product}")

    ann.run("CHR1:160", "GGT")
    check("multi-base ALT is declared, not translated",
          ann.mutation_type == "indel, not evaluated" and ann.snp_aa == "n/a",
          f"type={ann.mutation_type} snp_aa={ann.snp_aa}")
    ann.run("CHR1:160", ".")
    check("ALT '.' is not translated",
          ann.snp_aa == "SNP nt not provided", f"snp_aa={ann.snp_aa}")
    ann.run("CHR1:160", None)
    check("ALT None still reports ref_aa",
          ann.ref_aa != "n/a" and ann.snp_aa == "SNP nt not provided",
          f"ref_aa={ann.ref_aa} snp_aa={ann.snp_aa}")

    # start codon handling: build a CDS whose start codon is GTG (a table 11 initiator)
    g2 = "GTG" + "AAA" * 9 + "TAA" + "T" * 60
    rec2 = SeqRecord(Seq(g2), id="CHR2", name="CHR2")
    rec2.features = [SeqFeature(FeatureLocation(0, 33, strand=1), type="CDS",
                                qualifiers={"gene": ["altStart"], "product": ["p"],
                                            "transl_table": ["11"]})]
    a2 = Annotation(gbk_list=[], validate=False)
    a2.gbk_dict = {"CHR2": rec2}
    a2.run("CHR2:1", "T")          # GTG -> TTG, still a valid initiator
    check("GTG start reads M; GTG>TTG stays M and is silent",
          a2.ref_aa == "M" and a2.snp_aa == "M" and a2.mutation_type == "silent mutation",
          f"{a2.ref_aa}1{a2.snp_aa} type={a2.mutation_type}")
    a2.run("CHR2:1", "A")          # GTG -> ATG, also a valid initiator
    check("GTG>ATG stays M and is silent",
          a2.snp_aa == "M" and a2.mutation_type == "silent mutation",
          f"{a2.ref_aa}1{a2.snp_aa} type={a2.mutation_type}")
    a2.run("CHR2:2", "A")          # GTG -> GAG, not an initiator
    check("start codon destroyed reports 'start codon lost'",
          a2.mutation_type == "start codon lost",
          f"{a2.ref_aa}1{a2.snp_aa} type={a2.mutation_type}")

    # overlapping CDS must be reported, not silently resolved
    g3 = "ATG" + "AAA" * 40 + "TAA"
    rec3 = SeqRecord(Seq(g3), id="CHR3", name="CHR3")
    rec3.features = [
        SeqFeature(FeatureLocation(0, 60, strand=1), type="CDS",
                   qualifiers={"gene": ["geneA"], "product": ["a"]}),
        SeqFeature(FeatureLocation(30, 90, strand=1), type="CDS",
                   qualifiers={"gene": ["geneB"], "product": ["b"]}),
    ]
    a3 = Annotation(gbk_list=[], validate=False)
    a3.gbk_dict = {"CHR3": rec3}
    a3.run("CHR3:40", "T")
    check("two overlapping CDS are both reported",
          a3.overlapping_features == 2 and "geneA" in a3.overlapping_genes
          and "geneB" in a3.overlapping_genes and a3.cds_overlap == "yes",
          f"n={a3.overlapping_features} genes={a3.overlapping_genes} flag={a3.cds_overlap}")
    a3.run("CHR3:10", "T")
    check("single-CDS position reports no overlap",
          a3.overlapping_features == 1 and a3.cds_overlap == "n/a",
          f"n={a3.overlapping_features} flag={a3.cds_overlap}")


# --------------------------------------------------------------------- section 4

def gbk(name, pattern="*.gbk"):
    hits = sorted(glob.glob(os.path.join(DEP, name, pattern)))
    return hits or None


def section4_controls():
    print("\n[4] canonical drug-resistance controls, H37Rv NC_000962.3")
    g = gbk("Mycobacterium_H37")
    if not g:
        return skip("H37Rv controls", f"{DEP}/Mycobacterium_H37 not found")
    ann = Annotation(gbk_list=g, validate=False)
    chrom = "NC_000962.3"
    rec = ann.gbk_dict[chrom]

    def gene_feature(name):
        return next(f for f in rec.features
                    if f.type == "CDS" and f.qualifiers.get("gene", [""])[0] == name)

    # katG S315T, Rv1908c, minus strand
    katG = gene_feature("katG")
    e = int(katG.location.end)
    p0 = (e - 1) - (3 * 314 + 1)
    ann.run(f"{chrom}:{p0 + 1}", "G")
    check("katG S315T (minus strand)",
          (ann.ref_aa, ann.aa_residue_pos, ann.snp_aa) == ("S", 315, "T")
          and ann.reference_base_code == "AGC" and ann.snp_base_code == "ACC"
          and ann.mutation_type == "nonsynonymous",
          f"got {ann.ref_aa}{ann.aa_residue_pos}{ann.snp_aa} {ann.reference_base_code}>{ann.snp_base_code}")

    # rpoB S450L, Rv0667, plus strand
    rpoB = gene_feature("rpoB")
    p0 = int(rpoB.location.start) + 3 * 449 + 1
    ann.run(f"{chrom}:{p0 + 1}", "T")
    check("rpoB S450L (plus strand, must not regress)",
          (ann.ref_aa, ann.aa_residue_pos, ann.snp_aa) == ("S", 450, "L"),
          f"got {ann.ref_aa}{ann.aa_residue_pos}{ann.snp_aa}")

    # alternative start codons must read M
    for name in ("dnaA", "recF", "gyrB"):
        try:
            f = gene_feature(name)
        except StopIteration:
            continue
        s = str(f.extract(rec.seq))
        if s[:3] == "ATG":
            continue
        p0 = int(f.location.end) - 1 if f.location.strand == -1 else int(f.location.start)
        ann.run(f"{chrom}:{p0 + 1}", "A")
        check(f"{name} start codon {s[:3]} reads M at residue 1",
              ann.ref_aa == "M" and ann.aa_residue_pos == 1,
              f"ref_aa={ann.ref_aa} residue={ann.aa_residue_pos}")


def section5_translation_oracle():
    print("\n[5] /translation oracle: ref_aa must equal the reference protein")
    targets = [("Mycobacterium_H37", None), ("Mycobacterium_AF2122", None),
               ("Brucella_ovis", None), ("Mycobacterium_orygis", None),
               ("NC_045512_wuhan-hu-1", None)]
    random.seed(101)
    for name, _ in targets:
        g = gbk(name)
        if not g:
            skip(f"{name} oracle", "reference not found")
            continue
        ann = Annotation(gbk_list=g, validate=False)
        n = ok = 0
        first = None
        for chrom, rec in ann.gbk_dict.items():
            cds = [f for f in rec.features
                   if f.type == "CDS" and 'translation' in f.qualifiers
                   and 'pseudo' not in f.qualifiers]
            if not cds:
                continue
            for f in random.sample(cds, min(150, len(cds))):
                tr = f.qualifiers['translation'][0]
                order = [int(p) for p in f.location]
                off = random.randrange(0, (len(order) // 3) * 3)
                p0 = order[off]
                ci = off // 3
                if ci >= len(tr):
                    continue
                ann.run(f"{chrom}:{p0 + 1}", "A")
                if not ann.is_coding or ann.overlapping_features > 1:
                    continue          # ambiguous position, excluded by design
                n += 1
                want = "_" if tr[ci] == "*" else tr[ci]
                if ann.ref_aa == want and ann.aa_residue_pos == ci + 1:
                    ok += 1
                elif first is None:
                    first = (f"{chrom}:{p0+1} {f.qualifiers.get('gene',['?'])[0]} "
                             f"got {ann.ref_aa}{ann.aa_residue_pos} want {want}{ci+1}")
        if n == 0:
            skip(f"{name} oracle", "no checkable CDS")
        else:
            check(f"{name}  {ok}/{n} match /translation", ok == n, first or "")


def section6_sample_vcfs():
    print("\n[6] sample VCF fixtures")
    fixtures = [
        ("Mallard 26G08434-002", "26G08434-002-original_zc.vcf",
         [os.path.join(VCFDIR, "modified_NL_not_annotated", "Modified_NL_isolate.gbk")], 0),
        ("Owl 26G08579-001", "26G08579-001-original_zc.vcf",
         [os.path.join(VCFDIR, "owl_25-003495-001", "owl_25-003495-001.gbk")], 0),
        ("MTBC0 26G09227-017", "26G09227-017_zc.vcf", gbk("mtbc0_v1.1") or [], None),
        ("AF2122 26G08640-001", "26G08640-001_zc.vcf", gbk("Mycobacterium_AF2122") or [], None),
    ]
    for label, vcf_name, gbks, expect_minus in fixtures:
        vcf = os.path.join(VCFDIR, vcf_name)
        if not os.path.exists(vcf) or not gbks or not all(os.path.exists(g) for g in gbks):
            skip(label, "vcf or gbk not found")
            continue
        ann = Annotation(gbk_list=gbks, validate=False)
        total = coding = minus = bad = 0
        first = None
        for line in open(vcf):
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c[3]) != 1 or len(c[4]) != 1 or c[4] not in "ATGC":
                continue
            if c[0] not in ann.gbk_dict:
                continue
            total += 1
            ann.run(f"{c[0]}:{c[1]}", c[4])
            if not ann.is_coding:
                continue
            coding += 1
            if ann.direction == "reverse gene":
                minus += 1
            # invariant: reported codon must translate to the reported reference aa
            expect = ann._translate_codon(ann.reference_base_code,
                                          "11", ann.start_codon)
            if ann.ref_aa not in (expect, "unfound_ref_AA"):
                bad += 1
                if first is None:
                    first = f"{c[0]}:{c[1]} codon {ann.reference_base_code} ref_aa {ann.ref_aa} vs {expect}"
        ok = check(f"{label}  {coding}/{total} coding, codon/ref_aa consistent",
                   bad == 0, first or f"{bad} inconsistent")
        if expect_minus is not None:
            check(f"{label}  minus-strand CDS count == {expect_minus}", minus == expect_minus,
                  f"got {minus}")
        elif ok:
            print(f"        ({minus} minus-strand coding SNPs — these are the ones 3.35 got wrong)")


def main():
    print("=" * 78)
    print(f"vsnp3_annotation regression tests   (module version "
          f"{__import__('vsnp3_annotation').__version__})")
    print("=" * 78)
    section1_synthetic()
    section2_spliced()
    section3_boundaries_and_noncoding()
    section4_controls()
    section5_translation_oracle()
    section6_sample_vcfs()
    print("\n" + "=" * 78)
    print(f"PASS {len(PASS)}   FAIL {len(FAIL)}   SKIP {len(SKIP)}")
    if FAIL:
        print("failed:")
        for f in FAIL:
            print(f"  - {f}")
    print("=" * 78)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
