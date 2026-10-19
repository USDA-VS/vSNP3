# vSNP3 3.37 (unreleased) — changes since 3.36

Version 3.37. The literal was bumped with `internal/bump_version.sh 3.37`
(`bin/vsnp3_version.py`, `conda_build/meta.yaml`, `README.md`,
`container/Singularity.def`, every `internal/*.sh` wrapper) plus
`container/README.md` by hand; `tests/test_release_consistency.py` checks they
agree. The `meta.yaml` sha256 still has to be regenerated against the 3.37 tag.

Numbered items refer to the October 2026 notes. Items 1–25 were validated
against 3.36 on 2026-10-15 (2026-08-05 HPC validation dataset); where that
validation found an item inaccurate, the item below has been corrected and the
fix is described under "Corrections after the 2026-10-15 validation" (items
26–37, deviations D-01 to D-12 of the validation report).

## Samples are never dropped without being named (notes 6 and 8)

1. **A sample with no SNP passing the thresholds is now analysed, not dropped.**
   `vsnp3_group_on_defining_snps.py` used to append such a sample to
   `vcf_bad_list` as "Dataframe Empty" and leave it out of every group, table and
   tree. The step 2 summary then counted it under "VCF files in this run were
   corrupt" while the validation table said 0 invalid, and the sample was nowhere in
   the groupings table. That is the "3 corrupt / 0 failed / 413 of 416 listed"
   screenshot. B26-0290 (every call AC=1) is this case. Such a sample now sits in
   the All group and under "Group Not Found", normalized against the group's
   positions like any other sample: reference call where it matches, IUPAC code
   where it is mixed, `-` where it has no coverage. The summary lists each one with
   counts (records, zero-coverage, called, mixed AC=1, indel, single-base AC=2 and
   why those failed) and flags the mixed-sample signature: no homozygous SNP at all
   and nearly every call AC=1.
2. **A sample identical to the reference is analysed.** A header-only VCF was
   excluded as "empty" and a VCF holding only zero-coverage records was excluded as
   "no MQ field in INFO on any record". Both are now valid samples whose row is the
   reference call everywhere they have coverage. The MQ guard still applies when a
   file has called records and none carries MQ. Step 1 can now write such a VCF
   when freebayes calls no variant at all (item 36).
3. **Sample accounting section in the step 2 summary.** Every VCF found is in
   exactly one row: analysed (with the no-SNP subset), or not analysed under
   "failed validation", "could not be parsed", "removed by remove_from_analysis" or
   "duplicate sample name", each with its reason. The old "N VCF files were corrupt"
   line is gone; if the rows ever fail to add up the report says so.
4. **Step 1 says it too.** When no SNP passes the step 2 thresholds (single-base,
   AC=2, QUAL > 150, MQ ≥ 56) the step 1 report carries a notice: mixed calls
   present (likely mixed or contaminated sample), AC=2 calls that all fail QUAL or
   MQ, or no call at all (matches the reference where covered). As first written
   the notice was decided on "Quality SNPs", which counts any AC, so it could not
   appear for a typical mixed sample; see item 35.

## Reference type without a metadata file (note 1)

5. `vsnp3_step2.py -t <type>` on a reference directory with no `*meta*.xlsx`
   died with `UnboundLocalError: metadata_df`. The worksheet was only read inside
   the defining-SNP branch of `Group.__init__`; the same crash hit every
   `-abs_pos/-group` run. Metadata is now read on its own, is optional, and an
   unreadable or missing file is reported in the terminal and the summary while the
   run continues with VCF file names as sample names.
6. `vsnp3_reference_options.py` tested the whole path against `meta` and `remove`,
   so a reference directory anywhere under a folder with "meta" in its name (a
   "metagenomics" project, or my test directory) silently lost its defining-SNP
   file. The tests now run on the file name.

## Non-ASCII characters (note 7)

7. Every vSNP3 module forces `LC_ALL=C`; under it Python's text-mode `open()`
   decodes as ASCII on Linux, so an en dash or accented letter in a `##reference`
   path raised `UnicodeDecodeError`, the validator called the file corrupt, and the
   "auto-fix" rewrote the user's input in place. VCFs are now read as bytes, non-
   ASCII bytes are dropped from the in-memory copy and counted, the input file is
   never touched, and the summary carries one line ("N files contained non-ASCII
   bytes, ignored when read") with the names in the validation log. There is no
   longer a per-file "Fix: Removed ... character(s)" list.

## Run time reporting (note 5)

8. The summary now shows **Started**, **Finished** and **Total run time**, computed
   as the last thing before the file is written, plus a table of time per phase
   (validate, parse, archive, sort defining SNPs, select positions and annotate,
   group selection, normalize groups, trees and tables, report). 3.36 computed the
   run time, then ran `conda list <pkg>` twice for each of 13 packages (26 conda
   launches) and the version checks before closing the file. One `conda list
   --json` call replaces them, with `importlib.metadata` as the fallback.
9. A crash inside the HTML-tree child process hung the whole run: the wrapper only
   put a result on the queue on success and the parent waited on it without a
   timeout. The child now reports success or failure, a failure is logged per group
   and the run continues, and the 15-minute timeout still applies. A failure
   inside `html_tree` itself was still only printed (it returned 1, which counted as
   success) until item 33; failures are now listed at the end of the tree phase
   and in the step 2 summary.

## Speed (note 3)

10. **Validation** reads only the header and the first 20 records of each file,
    in a thread pool, and checks the reference from the `##contig` header lines
    instead of the first record's CHROM in a random sample of 50 files. It is
    deterministic. A minority of files on another reference are excluded and
    listed instead of aborting the run. The files kept must be a strict majority
    of the readable files; otherwise, an even split included, the run aborts as a
    mixed set (as first written an even split did not abort; see item 29).
11. **Reading** is one pass per file. 3.36 read each VCF three or four times: a
    header check, a line-by-line rewrite to `vcf_normalized/` with a dozen regex
    substitutions per line, a pandas read of that copy, a second rewrite for header
    notes, a copy to `vcf_starting_files/`, and a copy into the zip. It then split
    the INFO column into a dict three times per record for AC, DP and MQ. Now:
    binary read, normalization in memory with the regex work done only on lines
    that contain quotes or spaces, pandas straight from memory reading the eight
    fixed columns (the sample column was carried through every later stage and
    never used), and AC/DP/MQ by one vectorized regex each. The normalized copy is
    written only for `-f/--fix_vcfs` and `--debug`; `vcf_starting_files.zip` is
    built directly from the source files; `dictionary_of_dataframes.pickle` is
    written only with `--debug` (3.36 wrote the whole data set to disk on every run
    and deleted it a moment later). Measured on the six AF2122 test VCFs (113,331
    records) on a local SSD, single process, best of three: 3.36 reader 0.80 s,
    3.37 reader 0.22 s, a 3.6x speedup with identical frames. On a network file
    system the removed extra reads, writes and deletes of every input file should
    matter more than that ratio suggests; the per-phase table in the summary will
    show where the remaining time goes on a large run. The process pool is capped
    at the number of files and prints progress.
12. The heterozygous IUPAC mapping in `make_groupings` is vectorized (was a
    per-row `iterrows` loop with a try/except).
13. *Corrected after validation.* A record with more columns than the header is
    read normally: only the eight fixed columns are read, so the extra fields are
    ignored, and nothing is skipped or counted (3.36 also read it, truncated to ten
    fields). As first written this item said such a record was "skipped and
    counted"; it was not, and when the record was the first data line pandas
    rejected the whole file (see item 28). A line that is not a record at all — no
    numeric POS, such as a stray `#` line, a second VCF header or garbage — is
    dropped and counted under "Records that could not be read and were skipped"
    instead of becoming a record at position 0, as 3.36 parsed it. Blank lines are
    not counted (item 30).

14. **The tree-building phase no longer scales as the cube of the sample
    count.** The step 2 phase table showed "FASTAs, RAxML trees, tables and HTML
    trees" at 4 min 13 s of a 5 min 12 s run on 421 Brucella suis samples (81%),
    and the Mallard runs spend hours there. Measured on one synthetic group at a
    time (`tests/regression/bench_tree_phase.py`, below), this Mac, sequential
    `raxmlHPC`, best of one:

    | group | RAxML (unchanged) | tables, before → after | HTML tree, before → after |
    |---|---|---|---|
    | 100 samples × 2,000 positions | 28 s | 4.9 s → 2.1 s | 111 s → 0.9 s |
    | 200 × 3,000 | 101 s | 12.8 s → 5.5 s | **960 s** → 2.5 s |
    | 400 × 4,000 | 343 s | 33.8 s → 12.9 s | not run (≈ 2 h by the cube) → 5.4 s |
    | 1,000 × 8,000, HTML tree alone | | | → 34 s |
    | 2,000 × 8,000, HTML tree alone | | | → 124 s |

    960 s is past the 900 s timeout: in 3.36 and until now, a group of about 190
    samples or more with `-html_tree` was reported as "HTML tree generation timed
    out after 15 minutes" and got no HTML tree at all. The 421-sample all-VCF
    group would have been one of them. Every output is identical, see the data
    check at the end of this section. The changes, largest first:
15. **`vsnp3_html_tree.py` builds each tree from a bit matrix.** The tree is the
    same bottom-up agglomeration -- merge the pair sharing the most SNP
    descriptors, first pair in list order winning a tie -- but the pairwise shared
    counts come from one matrix product instead of a Python set intersection per
    pair per round, and a merge updates only the new node's row and the rows
    whose best partner it consumed. Descriptors with the same sample pattern are
    one weighted column. A bootstrap replicate is the same matrix with the drawn
    positions as weights (the draw itself is unchanged, same generator and seed),
    so none of the 100 replicates builds a resampled DataFrame. Clades are keyed
    by a bitmask of leaf names instead of a `'|'`-joined string rebuilt by a full
    leaf walk at every node. 102 tree builds per group went from 104 s to 0.9 s at
    100 samples. Verified node by node against the old builder (structure, names,
    SNP lists, Excel row indexes, identical-sample groups, bootstrap counts) on
    400 random tables with ties, no-calls, N, IUPAC codes, identical samples and
    numeric names, as well as on the AF2122 set and the synthetic groups above.
16. **Plotly annotations and traces are added in one call per figure.**
    `add_annotation()` re-validates and deep-copies every annotation already on
    the figure, so one label per sample was quadratic: 35 s at 100 samples. The
    figure holds the same traces and annotations in the same order; the HTML is
    byte-identical apart from plotly's per-render uuid.
17. **The HTML tree is built from the frame step 2 just wrote, not from the
    workbook read back through openpyxl** (0.8 s at 100 samples, 6.4 s at 400,
    longer than the whole tree now takes). `Tables` keeps the cascade1 frame as
    `read_excel` would return it (`as_read_back`: index as the first column, the
    default NA strings as missing, the coverage column as `excel_formatter`
    rewrites it) and `html_tree()` takes it as `df=`; the benchmark checks the two
    frames agree, and the CLI `vsnp3_html_tree.py -f table.xlsx` still reads the
    file and produces the same HTML.
18. **The per-column table statistics are one numpy pass.** The three cascade
    sorts counted SNPs per column, the first SNP row and the longest vertical run
    in nested Python loops that re-indexed the root cell for every element
    (463,000 `iloc` calls for 200,000 cells): 7.0 s → 0.1 s at 400 × 4,000.
19. **Workbooks are written cell by cell through xlsxwriter directly** rather
    than `to_excel`, which built a cell object per value before calling the same
    writer; the JSON round trip each table makes on the way (it rounds floats to
    10 digits and must stay) is done in memory instead of through `df.json` in
    the output directory. Three workbooks: 26.5 s → 12.2 s at 400 × 4,000. The
    xlsx bytes differ (shared-string order), the cells do not. As first written,
    column A (sample names and the root, MQ and annotation labels) lost the format
    `to_excel` gave the index; it is restored in item 27.
20. **RAxML runs are no longer oversubscribed.** Each run keeps `-T 4` (same
    binary, same command, same tree), but a semaphore allows only cores ÷ 4 runs
    at once, counting 1 thread for a non-PTHREADS build, which ignores `-T`. A
    16-core machine used to run 13 × 4 = 52 busy-waiting threads. The pool is not
    shrunk, so a group whose RAxML has finished goes on to its tables while others
    wait. `VSNP3_RAXML_THREADS` overrides the 4 and `VSNP3_MAX_CPUS` caps the cores
    counted. Groups are submitted largest alignment first: the all-VCF group is
    defined last and is the biggest, so the phase used to end with it running
    alone. The RAxML version banner is read once per process instead of once per
    group.
21. Every walk of the HTML tree is recursive, and a chain-like group is as deep
    as it has samples, so the default recursion limit failed such groups past a
    few hundred samples with "Error: maximum recursion depth exceeded" (printed,
    not reported). The limit now scales with the sample count.

    Data check on all of the above: `tests/regression/compare.py` on the AF2122
    set with `-a -html_tree`, before against after: **data identical**, the four
    `_position_tree.html` files byte-identical after uuid normalization, and two
    consecutive runs byte-identical. `tests/regression/bench_tree_phase.py` is
    new: it generates a synthetic group of any size and times each sub-step of
    this phase against any `bin/` directory, caching the RAxML tree so two
    versions are timed on the same tree.

Data check: `tests/regression/compare.py` on the AF2122 test set, 3.36 against
3.37 with `-a -html_tree`: **data identical**.

## Step 1 usability verdicts (note 4)

22. "FASTQ Usability" and "Reference Usability" are no longer written to the stats
    workbook (`HTML_ONLY_KEYS` in `vsnp3_file_setup.py`), so merged Excel files
    carry measurements only.
23. In the HTML report each verdict tile shows the measurement that set it, e.g.
    "Questionable — average depth 35.2x < 40.0x". The limits now live in
    `vsnp3_version.py` (`FASTQ_USABILITY`, `REFERENCE_USABILITY`) and were loosened:

    | | 3.36 Poor | 3.36 Questionable | 3.37 Poor | 3.37 Questionable |
    |---|---|---|---|---|
    | average depth below | 40x | 70x | 20x | 40x |
    | duplicate reads above | 80% | 10% | 80% | 30% |
    | reads passing Q20 / Q30 below | 50% Q20 | 70% Q30 | 50% Q20 | 70% Q30 |
    | reference covered below | 95% | 97.5% | 90% | 97% |
    | zero-coverage above | 1% | 0.5% | 10% | 3% |
    | quality SNPs, % of reference, above | 0.09% | 0.03% | 0.09% | 0.03% |

    The 3.36 "zero-coverage above 1%" rule meant any sample below 99% coverage was
    Poor, whatever the coverage rule said. Adjust the dicts if these are still not
    right; nothing else needs to change.

## Documentation (note 2)

24. README (Step 2 section and More reading) and
    `docs/instructions/additional_tools.md` link the
    [Phylogenetic Tree HTML Viewer](https://usda-vs.github.io/phylogenetic-tree-html-viewer/)
    for viewing and sharing vSNP trees.

## Tests

25. `tests/test_step2_inputs.py`: 51 checks on synthetic VCFs covering every item
    above that can be tested without RAxML (validation categories and determinism,
    inputs left unmodified, reader output and statistics, no-SNP sample handling in
    `make_groupings`, the no-metadata path, the child-process wrapper, and the
    workbook/report split). `tests/test_smoke.py` still passes (113 checks).

## Corrections after the 2026-10-15 validation

The validation of 3.37 against 3.36 found every SNP call, table cell and tree
identical on the validation dataset, and listed deviations D-01 to D-12 (report
section 7). Each fix below has a check in `tests/test_step2_inputs.py` that fails
on the code as validated, except items 26 and 36, checked by
`tests/test_release_consistency.py` and by a step 1 run respectively.

26. **D-01, D-07: version and release folder.** The version literal said 3.36, so
    every 3.37 output, log and report called itself 3.36. Bumped to 3.37 as
    described at the top. The release folder gained `LICENSE` (declared in
    `meta.yaml`; GPL-3.0, from the GitHub repository), `.gitignore` (from the
    development repository, with the local working-notes entry the release check
    requires) and
    `internal/bump_version.sh` (the release script `bin/vsnp3_version.py` refers
    to). `tests/test_release_consistency.py` now also works in a folder that is not
    a git checkout: it scans every file `.gitignore` lets through (it used to see
    no files at all there and pass without reading anything) and reads
    `.gitignore` itself; commit messages and authors can only be checked in the
    repository, and it warns that it did not. Six wrappers that carried 3.19–3.35
    now say 3.37, as the bump script sets them.
27. **D-02: column A of every workbook keeps its 3.36 format.** Sample names and
    the root, MQ and annotation labels are bold, thin-bordered, centred and
    top-aligned, as `to_excel` wrote the index. The annotation label no longer
    takes the annotation row's rotated blue format. Cell values and every other
    format are unchanged (`Tables.excel_formatter`, `Tables._write_cell`).
28. **D-03: an extra field on the first data line no longer rejects the file.**
    3.37 as validated excluded such a file under "could not be parsed" ("Number
    of passed names did not match number of header fields"), where 3.36 analysed
    it with correct calls. The reader now names only the eight columns it reads,
    so the line is read like any other and its extra fields are ignored; item 13
    is corrected to say so. The 3.37 test that seemed to cover this passed only
    because its line had no numeric POS.
29. **D-04: an even split between references aborts.** Three B. suis and three
    M. bovis files kept the half whose contig names sort higher. The files kept
    must now be a strict majority of the readable files (`kept * 2 > readable`);
    otherwise step 2 stops with "Only N of M readable files are on the most
    common reference ..., which is not a majority" and nothing is analysed. A
    minority on other references is still excluded and named as before.
30. **D-05: blank lines are not reported as skipped records.** A VCF ending in a
    blank line was analysed with identical calls but listed under "Records that
    could not be read and were skipped". Blank lines are no longer counted as
    records.
31. **D-08: `VSNP3_BOOTSTRAP`, an optional RAxML bootstrap mode, documented.**
    This was in 3.37 but in no item above. It comes from vsnp_gui's
    `bootstrapfix.py` patch (the `KL_BOOTSTRAP` block in
    `bin/vsnp3_fasta_to_snps_table.py`; vsnp_gui sets the variable when its Step 2
    bootstrap field is above 0). With the environment variable set to N > 0,
    every RAxML tree of the run (each group, and the all-VCF tree with `-a`) is
    built by RAxML's rapid bootstrap analysis, `-f a -x 7777 -N N` (seed
    `-p 456123` as before), and the tree kept is `RAxML_bipartitions`, the
    best-scoring tree of that analysis with bootstrap support on its internal
    nodes. Unset, empty, 0, negative or not a whole number: best tree only, exactly
    as 3.36 (newick byte-identical on six test sets). This is a different analysis,
    not the default tree with labels added: on a 40-sample test set the topology
    differed from the default run (Robinson–Foulds 4 of 38 splits) and so did the
    table row order. The tree stage takes much longer. **The 2026-10-15
    validation covered the default only: leave `VSNP3_BOOTSTRAP` unset for
    validated runs until this mode is validated in its own right.** The value is
    read by `bootstrap_replicates()`; the misleading "same ML topology" comment is
    corrected.
32. **D-08: the step 2 summary states the tree mode.** Under the thresholds:
    "RAxML trees: best tree only (VSNP3_BOOTSTRAP not set)", or the value set and
    how it was read, or a highlighted box "RAxML trees: rapid bootstrap, N
    replicates (VSNP3_BOOTSTRAP=N)". The line is new in every summary.
33. **D-10: a failed HTML tree is recorded, not only printed.** `html_tree`
    catches its own errors, prints them and returns 1, and the child-process
    wrapper took the 1 as a normal result, so `html_tree_failures` stayed empty.
    Step 2 now calls it with `raise_errors=True`, so the exception and its message
    reach the wrapper, and any non-zero return is a failure as well. Failures are
    listed at the end of the tree phase and in a "HTML tree not written" box in the
    summary (shown only when there is one); the group's tables and RAxML tree are
    unaffected. The HTML-tree step is now `Group.html_tree_build`. The command line
    `vsnp3_html_tree.py` still prints and exits 1.
34. **D-11: sample names such as 'NA', 'None' or 'null' stay names in the HTML
    tree.** The frame handed to `html_tree` mirrored `read_excel`, which reads
    those strings as missing: the leaf was labelled "nan" and two such samples
    shared one clade key (104% support). Both that frame (`Tables.as_read_back`)
    and the workbook read on the command line (new `vsnp3_html_tree.read_table`)
    now read the sample-name column with `keep_default_na=False`; every other
    column is read exactly as before. Tables were already correct.
35. **D-12: the step 1 no-SNP notice uses the step 2 definition.** "Quality SNPs"
    counts single-base calls with QUAL > 150 of any AC, so a mixed sample whose AC=1
    calls are of high quality (the B26-0290 case; CT14b had "Quality SNPs 20") never
    got item 4's notice. Step 1 now counts what step 2 will take as SNP positions —
    single-base, AC=2, QUAL > 150, MQ ≥ 56, QUAL and MQ truncated, MQM read as MQ,
    as `Group._passing` — and gives the notice when that count is 0
    (`step2_snp_counts` and `no_snp_notice` in `bin/vsnp3_zero_coverage.py`; a test
    checks the count against `Group._passing`, boundaries included). The terminal
    shows the count. "Quality SNPs" itself, its tile and the reference-usability
    limit it feeds are unchanged. On the validation dataset every sample has at
    least 414 (M. bovis) or 69 (B. suis) such SNPs, so none gains a notice;
    B26-0290 and CT14b have 0 and now get the mixed-sample notice.
36. **D-09: step 1 no longer aborts on a sample identical to the reference.**
    Reads with no difference from the reference gave a freebayes VCF with a header
    and no record, and step 1 stopped with "freebayes output contains no variant
    records", no VCF written and the intermediate files left behind. A VCF with its
    `#CHROM` header and no record is now accepted (the reads are already required
    to have mapped); step 1 goes on to a `_zc.vcf` of zero-coverage records only,
    notes it in the run log, and the report carries the "matches the reference"
    notice. A file without a header is still a failed run. The nanopore path
    (bcftools) is unchanged and still stops in that case.
37. **Tests.** `tests/test_step2_inputs.py` grows from 51 to 78 checks: sections
    [8]–[12] (reference majority, workbook column A, NA-like names, HTML-tree
    failures, step 1 notice) and new reader cases (extra field on the first and a
    later line, blank lines). `tests/regression/bench_tree_phase.py` compares the
    handed-over frame with `read_table` where the code under test has it.

Not changed: D-06 (a truncated VCF is analysed as if complete; 3.36 behaves the
same) remains open for a later release. Until then, checksum VCFs copied between
systems.

## Corrections after the 2026-10-15 accuracy study

The accuracy study (record VAL-vSNP3-3.37-20261015-A2) scored 3.37 against known
truth and documented three limitations: L-1 (SNPs in repeated sequence are not
called; by design, unchanged), L-2 and L-3 below, and deviation D-06 was still
open.  The version stays 3.37.

38. **L-2: a SNP inside a multi-base record is read from that record.** freebayes
    sometimes reports a single-base change inside a longer record that starts at
    an earlier position, usually because a nearby read carries an indel: for
    example REF=CAC ALT=CCC at NC_002945.4:2324531 for an A>C change at 2324532
    in 13-1941.  Step 2 looked records up by their start position only, found
    none at the SNP's own position and wrote the reference base, so the sample
    looked one SNP different from relatives that share the change.  The reader
    now keeps the CIGAR of multi-base records, and `make_groupings` reads the
    sample's base at every group position inside such a record from it (M, X
    and = align a base, D deletes one, I inserts).  The call then follows the
    usual rules: the base at QUAL >= 150, an IUPAC code when the record is mixed
    (AC=1), N at QUAL 50-150, the reference below 50.  It is N where the base is
    deleted, where the record has no usable CIGAR and REF and ALT differ in
    length, or where the record's alleles disagree.  A base the record leaves as
    the reference stays the reference.  Positions chosen for the tables and the
    first position of each record are unchanged.  A position can leave a group's
    table as a result: with the base read correctly every sample may carry the
    same base.  The summary report lists the cells read this way and the
    positions that left a table.  On the 2026-08-05 validation samples this
    affects one position (13-1941 at 2324532, which now reads C and, carried by
    every M. bovis sample, leaves the four tables).  On 155 M. tuberculosis
    samples 826 table cells change, at 383 sample-positions in 150 samples: 55
    now show the sample's other base, 6 an IUPAC code and 765 N (763 of them
    where the sample has a deletion); 12 group positions leave a table.  3.36 and
    earlier showed the reference base at all of them.  The accuracy study's
    estimate of 499 sample-positions compared REF and ALT base by base, which an
    indel in the record misaligns; read through the CIGAR, most of those
    positions hold the reference base.
39. **L-3: a mixture with passing SNPs is flagged.** The 3.37 mixed-sample notice
    and summary flag fired only when no SNP passed the thresholds; a 60/40
    mixture in which one site happened to be called homozygous was not flagged.
    A sample is now also flagged, in the step 1 report and in a "Possible mixed
    samples" section of the step 2 summary, when its confident mixed calls
    (single-base AC=1, QUAL >= 150, MQ >= 56) number at least 10 and outnumber
    its passing SNPs (`vsnp3_version.mixed_signature`, `MIXED_MIN_AC1`,
    `MIXED_AC1_RATIO`).  The limits come from the validation data: known mixtures
    had 20 to 87 confident mixed calls per passing SNP, clonal samples at most
    0.41.  Flagged samples are analysed as before.
40. **D-06: a truncated VCF is excluded.** The validator reads only the start of
    each file, so a VCF cut off part-way (an interrupted copy, a full disk) was
    analysed as if complete, with every SNP after the cut read as the reference.
    A file that does not end with a newline ends part-way through a record; it is
    now excluded and listed under "Not analysed" as truncated.  One byte is read
    per file.  A cut that falls exactly on a line boundary cannot be seen this
    way.
41. **Tests.** `tests/test_step2_inputs.py` grows from 78 to 103 checks: sections
    [13] (multi-base records: every CIGAR case, mixed and low-quality records,
    deletions, insertions, records without a CIGAR, two alleles, and a position
    that leaves the table), [14] (the mixed-sample rule, step 1 and step 2
    counting the same calls) and [15] (truncated files, including the limit).

## Long-read (nanopore) data

Until now a nanopore sample was aligned with minimap2 and called with bcftools
mpileup/call at default settings, and 100 was added to every QUAL.  bcftools'
QUAL is the same quantity as freebayes' (the Phred-scaled posterior that the
site is a variant) but it cannot exceed about 228, because the genotype
likelihoods behind it are capped at PL 255; freebayes' grows with the evidence
into the thousands.  So a strong nanopore call sat at 225-228 (328 after the
offset), a moderate one anywhere below, nothing could reach the 300 consensus
threshold of `vsnp3_vcf_merge_to_fasta.py`, and the 50/150 bands of step 2 did
not mean for a nanopore VCF what they mean for an Illumina one.  The version
stays 3.37; the Illumina path is unchanged (items 42-48 below, validated in
record VAL-vSNP3-3.37-20261015-A4).

42. **A nanopore VCF now carries a QUAL on the freebayes scale.** The new
    `bin/vsnp3_longread.py` runs `bcftools mpileup` with its nanopore profile
    (`-X ont`: no BAQ, bases at Q5 or better, base quality capped at 30, no
    indel calling; written out as the same flags on a bcftools older than 1.17),
    `-d 10000` so a deep site is not subsampled, and `bcftools call -mv
    --ploidy 2`, so AC=1 and AC=2 mean what they mean for freebayes.  bcftools
    decides which sites are variants and their genotype.  QUAL is then
    replaced by the **evidence QUAL**: the Phred-scaled posterior probability
    that the site is reference-only, computed from the high-quality read
    counts bcftools used (INFO/DP4), a per-read error rate measured from the
    sample's own alignment (`samtools stats`: mismatches and indels per
    aligned base, clamped to 0.2%-25%) and the standard variant prior of 0.001.
    It is not capped: a clean site scores about 10·log10((1-e)/e) per
    supporting read -- 18 at a 1.5% error rate, 27 at 0.2%, which is what
    freebayes gives a Q30 Illumina read -- less about 33 for the prior.  A
    nanopore read carries less evidence than an Illumina read, so more
    nanopore reads are needed for the same QUAL: about 10 clean reads for 150
    at a 1.5% error rate, where freebayes needs about 6.  The thresholds keep
    their meaning on both platforms, and the +100 offset (`NANOPORE_QUAL_OFFSET`)
    is gone.  Records with caller QUAL <= 20 or that are not single-base SNPs are
    removed, as before.  QUAL is never below 50: bcftools reported the site as
    a variant, so however thin the evidence the sample is uncertain there and
    step 2 shows N, never the reference base (in the validation, a 10x sample
    showed the reference at true SNPs supported by four or five reads before
    this floor was added).  The genotype comes from the same model: a site is
    mixed (AC=1) only when the mixed hypothesis beats homozygous ALT under the
    measured error rate, which at 30x needs a minor allele of about a fifth of
    the reads.  bcftools trusts each base's own quality (capped at 30), so four
    reference-looking reads among thirty, the usual systematic nanopore error,
    made it call a clonal sample mixed at a true SNP, an IUPAC code in the
    table and a count towards the mixed-sample flag.  Where the model decides
    otherwise AC and GT are rewritten and bcftools' genotype is kept as
    INFO/VSNP_CALLER_GT.  bcftools' own QUAL is kept as INFO/VSNP_RAWQ, the
    evidence QUAL before the floor and the limits as INFO/VSNP_EVQ.
43. **Two nanopore error modes are placed in the N band.** Read counts cannot
    see them: ALT support on one strand only although both strands are
    sequenced (at least 4 ALT reads and 2 reads per strand), and a call in or
    beside a reference homopolymer run of 6 bases or more, where the
    indel-prone reads align least reliably.  Such a record is given a QUAL
    between 50 and 149 whatever its evidence, with the reason in FILTER
    (`vsnp_strand`, `vsnp_homopolymer`) and INFO/VSNP_FLAG, so step 2 reads the
    sample as N there: uncertain, never a false SNP or a false reference call.
    Runs of 6 or more touch 0.18% of M. bovis AF2122 and 0.41% of B. suis 1.
    The limits are constants in `vsnp3_version.py` (`LONGREAD_*`).
44. **The VCF header says how it was made**: `##vsnp3_platform=nanopore`,
    `##vsnp3_longread_caller=...` (the bcftools command) and
    `##vsnp3_longread_model=...` (the measured error rate and its source, the
    prior, the per-read weight, the strand and homopolymer limits), plus
    `##INFO` and `##FILTER` definitions for the new fields.  Step 2 reads the
    platform line and names the long-read samples in its summary, with each
    sample's model line, and in the terminal; the thresholds themselves are
    applied unchanged.  Nothing in step 2 treats a nanopore sample differently.
45. **Alignment and BAM.** `minimap2 -x map-ont --secondary=no`, read group
    `PL:ONT` (it said ILLUMINA).  Long reads are sorted and indexed, not
    duplicate-marked: nanopore libraries are not PCR-amplified, and `samtools
    markdup -r` on single-end reads compares start positions only, so at depth
    it discarded reads that merely began at the same base.  The BAM keeps the
    `_nodup` name every wrapper looks for; the read counts come from `samtools
    flagstat`, the duplicate fields of the statistics say "not marked (long
    reads)", and the run log says so.  The unmapped reads are written out for a
    nanopore sample as for any other.  A bcftools VCF with a header and no
    record (a sample identical to the reference) is accepted, as the freebayes
    one is since item 36; it used to stop step 1.
46. **Statistics and report.** Every step 1 run, both platforms, now records
    `SNP Caller` (freebayes or bcftools with its version and options).  A
    nanopore run adds `Long-read Error Rate`, `Long-read QUAL Model` and
    `Long-read Calls Limited to N` (how many, one-strand and homopolymer), shown
    under SNPs in the HTML report.  The run log carries the bcftools version
    and each stage of the calling.  Platform detection is unchanged: `-p`, or
    single-end reads averaging more than 701 bases.
47. **`vsnp3_longread.py` on its own** calls and models a BAM
    (`-r reference.fasta -b sample_nodup.bam`) or applies the model to a VCF
    bcftools call has already written (`--vcf raw.vcf`), with `--error_rate`
    to override the measurement; `-d` keeps the raw bcftools VCF.
48. **Tests.** `tests/test_step2_inputs.py` section [16], 103 to 146 checks:
    the model's shape (grows with supporting reads, falls with reference
    reads, confident at 15/15, below the N band at 3 of 30, the closed form on
    the freebayes scale, no overflow at depth), the homopolymer and strand
    rules, a bcftools-style VCF through `apply_model` (counts, header lines,
    every record's QUAL, FILTER and INFO, the AD fallback, FORMAT untouched)
    and then through step 2's reader, `Group._passing`, the step 1 counters and
    the group reporter.  `tests/test_smoke.py` 117 to 120 (the new CLI).

## Step 2 summary layout and vSNP GUI compatibility (2026-10-16)

49. **The "Not analysed" list is at the end of the summary, once.** The step 2
    summary listed every file left out twice: with its reason and a red X under
    "Sample accounting" near the top, and again as a bare list of names under
    'VCF files removed from dataset using "remove_from_analysis.xlsx"' after the
    Groupings table (which also named only the remove-list files, not those that
    failed validation or shared a sample name). The bare list is gone. The red-X
    list, grouped by reason, now sits where it was, under "Not analysed (N)", and
    the Sample accounting table at the top keeps the counts with a note that the
    files are listed at the end. `HTML_Summary` no longer takes `removed_samples`.
    Regression on the six AF2122 VCFs with a one-sample remove list: data
    identical; the summary differs only by this move.
50. **`VSNP3_MAX_CPUS` caps every worker pool.** Pools were sized at
    `cpu_count / 1.2` from the machine's core count, which is not the
    allocation: an 8-core job on a 64-core node, or a session on a 128-core
    shared server, started dozens of processes each holding its share of the VCF
    frames, and a 1,951-sample step 2 in a 32 GB session was killed with no
    traceback. The vSNP GUI and its suite launcher (bioinformatic_diagnostic_tools)
    pass the session's core budget as `VSNP3_MAX_CPUS`, which until now only the
    RAxML slot limit in `vsnp3_fasta_to_snps_table.py` read; the GUI's install
    step rewrote the pool-size expression in `vsnp3_step2.py` and
    `vsnp3_group_on_defining_snps.py` to add the cap. `vsnp3_file_setup.worker_count()`
    now sizes both pools: `cpu_count / 1.2`, capped by `VSNP3_MAX_CPUS` when it is a
    positive whole number and by 32 (`DEFAULT_MAX_WORKERS`) otherwise, never more
    than the number of tasks. The default cap is the one the GUI's patch has
    applied since May 2026, so a 3.37 install behaves as the patched installs
    did; on a 48-core HPC node it means 32 readers instead of 40 for a phase
    that took 6 minutes of a 5.5-hour run. RAxML concurrency is unchanged (its
    own slot limit). The `Tree.cpu_count` attribute, assigned and never read, is
    removed.
51. **The GUI's install patches stand down on 3.37.** vsnp_gui's
    `deploy/vsnp3-patches/apply.sh` rewrites an installed vsnp3 after every conda
    install. Each of its content fixes was checked against this release: the
    passQ comma fix (no `float(passQ20)` here; `_usability` converts with a
    try), the CPU cap (item 50: both files now name `VSNP3_MAX_CPUS`, its skip
    condition), the optional-metadata fix (3.37 assigns `metadata_df = None`
    itself; the comment now carries the marker phrase the script checks, so it
    does not insert a second assignment), `strandfix.py` (recognises the upstream
    minus-strand correction and exits 0) and `bootstrapfix.py` (sees
    `VSNP3_BOOTSTRAP` and exits 0). The v3.16-only `.patch` set is skipped on
    any other version. So a 3.37 install is left byte-identical by `apply.sh`.

## Robustness fixes from the 2026-10-16 audit, step 2 (2026-10-16)

An audit of every script, the packaging and the two vSNP GUI repositories listed
64 findings (kept with the working notes as `vsnp3_audit_findings_2026-10-16.md`,
IDs H1-H12, M1-M29, L1-L23, G1-G8). The confirmed ones that make vSNP3 more
robust without changing what a current run produces were fixed. Every fix below
was checked with the regression harness on the AF2122 set (`compare.py`: data
identical to the run before the change), the step 2 input tests and the smoke
test; the group normaliser was also run before and after on 400 seeded random
synthetic groups with identical output.

52. **Newick leaf order is read with Bio.Phylo (H1).** The tree's leaf order,
    which fixes the row order of every table, came from a chain of regex
    substitutions over the Newick text. The one meant to strip branch lengths,
    `[0-9].*\.[0-9].*\n`, matched anywhere in a line, so a sample name containing
    a digit, a dot and a digit was deleted (`22-12345.1_zc.vcf` vanished) or cut
    and glued to the next leaf (`NZ_CP0123.1.vcf` and `B_zc.vcf` became
    `NZ_CPB_zc.vcf`). The filter that followed kept only names present in the
    alignment, so the samples left every table and the HTML tree with no message.
    Metadata names were safe (`.` becomes `_`); file-name-derived names were not.
    `Tables.leaf_order` now reads the tree with `Bio.Phylo` in file order, and a
    sample in the alignment that is not in the tree is an error, as the reverse
    already was. Row order is unchanged for names the regex did not damage.
53. **No phantom taxon without a GenBank file (H2).** The HTML tree took every
    table row except `MQ` and `annotation` as a sample. Without a gbk the table
    carries an all-NaN row named `no annotations`, which became a leaf differing
    from the root at every position. `NON_SAMPLE_ROWS` in `vsnp3_html_tree.py`
    names all three rows. Checked with a reference set holding no gbk: four
    occurrences of the phantom in the tree before, none after, data identical.
54. **A blank metadata cell no longer names a sample "nan" (H3).** `astype(str)`
    turned a blank cell into the string `nan`, so one blank renamed its sample
    and two blanks collided and dropped a sample under "duplicate sample name:
    nan". Rows with a blank cell in either column are skipped (the VCF file name
    then names the sample), counted, and listed in the summary.
55. **Group failures reach the summary; names RAxML cannot carry are made safe;
    samples with no determined base are left out of a group and named (H4).**
    A group whose RAxML or tables failed was printed in the terminal only, and
    the Groupings table listed its members as analysed; the summary now carries a
    "No table or tree" box and marks the group in the Groupings table. Two causes
    were removed: a file-name-derived sample name containing a space, quote,
    bracket, colon, comma or semicolon failed RAxML for the whole group (such
    names get `_` in their place, reported in a "Sample names changed" box;
    metadata names were already cleaned); and a sample whose every call at the
    group's positions is `-` or N made RAxML exit ("sequences that consist
    entirely of undetermined values"), taking the group with it. Such a sample is
    left out of that group, and the parsimony test is repeated without it until
    every remaining sample has a determined base at some retained position (a
    sample determined only at positions the test removes would fail RAxML just
    the same); it is named in the
    terminal and under "Left out of a group's tree and table" in the summary, and
    stays in every other group it belongs to. Checked on 1,200 seeded random
    groups against the code before the change: identical output for every group
    with no such sample, and for the others exactly the placeable samples.
56. **The parsimony test is O(positions) in memory (H5).** It extended two lists
    with one `abs_pos+ALT` string per sample per position, one of them never read,
    and counted them: about 27 GB for an all-VCF group of 10,000 samples by
    20,000 positions, inside the "Normalize groups" phase that took 1 h 10 min in
    the 2026-10-15 production run. A position is uninformative when every sample
    carries the same single-base call, which is "no sample differs from the first",
    kept as a boolean mask. The old test only ever removed positions whose shared
    call was one character long (it recovered the position by dropping the last
    character), and that condition is kept explicitly, so the result is unchanged.
57. **The defining-SNP worksheet is checked once after reading (M1).** A position
    header that is not `chrom:position` text or a blank group name used to crash
    the run deep in the grouping stage (`AttributeError: 'int' object has no
    attribute 'split'`, `KeyError: nan`), after the input VCFs had been deleted.
    They are now a clear `DefiningSnpsError` naming the cell, and the run exits
    before anything is removed (item 66). A group name given by more than one
    defining position used to be overwritten, so the first position's samples left
    the group silently; the group is now the union of both, its filter columns are
    merged, and the summary says so. Filter entries that are not `chrom:position`
    text never matched a record; they are now listed, as they are still ignored.
58. **One selection pass per sample (M2).** `_passing` was evaluated up to three
    times per sample over every record, zero-coverage ones included, each with two
    `str.len()` passes: about 36 minutes of a 10,000-sample run. The reader writes
    a `single_base` column once per file, the selection is made once and the
    filtered selection and the counts reuse it.
59. **A table split into several workbooks is reported as having no HTML tree
    (M3)** instead of the tree being silently absent for exactly the largest
    groups.
60. **The HTML-tree timeout scales (M4).** It was a constant 900 s; at the
    measured growth (34 s for 1,000 samples, 124 s for 2,000) the all-VCF tree
    of a 10,000-sample run would always have timed out. `VSNP3_HTML_TREE_TIMEOUT`
    (seconds) sets it; otherwise it is 900 s or twice the fitted time for the
    group's size, whichever is larger.
61. **`count_leaves` is iterative (M5).** It recursed through `sum()` over a
    generator at every level, which the raised recursion limit does not protect:
    a chain of 3,000 samples still raised `RecursionError`.
62. **A remove-list entry that matches no file is reported (M6)**, in the terminal
    and in a summary box, instead of silently leaving the sample in every tree.
63. **The reference is detected by exact FASTA header ID (M7).** When no `-t` is
    given, the first CHROM of the first VCF was matched as a substring of the
    joined FASTA headers, so a short or generic contig name picked the first
    directory in sorted order whose header text contained it. A directory whose
    header ID equals the contig wins; the substring test remains only as the
    fallback, and either an ambiguous or a fallback choice is printed and shown in
    the summary.
64. **`vsnp3_excel_merge_defining_snps.py` (M8)** stripped the file name from the
    first `.` of the whole path, so `./old/a.xlsx` and `./new/a.xlsx` both became
    `_expanded.xlsx` and the merge compared a file with itself. It uses the
    basename, and a NaN cell or a non-text position cell is an error naming the
    column rather than a TypeError.
65. **HTML trees are half the size (M24).** Each file embedded plotly.js twice
    (once per figure) and also loaded `plotly-latest` from the CDN, a URL frozen at
    plotly.js 1.58 and deprecated, which needed network access and was then
    overwritten by the embedded bundle. The bundle is embedded once, with the
    first figure; the CDN tag is gone. 9.3 MB became 4.7 MB per group on the test
    set.
66. **Smaller step 2 items.** A sample name starting with `=` is written as text,
    not as an Excel formula (L6). The step 2 argparse defaults and the summary's
    threshold line come from `vsnp3_version.py` (`MQ: >=56`, as the selection
    tests), the density defaults are written once, and `vsnp3_vcf_merge_to_fasta.py`
    imports its thresholds (L7). "All_VCFs is available" is printed only when the
    all-samples group produced a table (L8). The two `except ValueError` clauses in
    the call rewrite print and record the error instead of hiding it behind
    `--debug`, and the summary lists the pairs (L9). The input VCFs are deleted at
    the end of the run, after the summary, not before the grouping stage, so a
    worksheet mistake no longer leaves a directory holding only the zip (L10).
    `vsnp3_vcf_merge_to_fasta.py`: `rstrip('_zc')` stripped characters, not the
    suffix, and a multi-contig reference is refused with a plain message (L11).
    `vsnp3_vcf_annotation.py` names its temporary files after the VCF and strips
    only the extension of the output name (L4). The "no VCF files" message names
    the pattern actually used (L13).
67. **UTF-8 reads on the step 2 side of H6.** GenBank files are parsed through an
    explicitly opened UTF-8 handle (an accented author name under `LC_ALL=C` was
    a `UnicodeDecodeError` on Linux, and the path form left a `ResourceWarning`),
    and `vsnp3_reference_options.py` and `vsnp3_path_adder.py` read and write the
    reference-paths file and FASTA headers as UTF-8. The step 1 side is item 78.
68. **For the vSNP GUI (G3, G5).** Step 2 writes
    `vSNP_step2_sample_accounting-<stamp>.json` beside the summary: every file
    found, analysed (with the no-SNP subset), not analysed with reason, remove-list
    entries that matched nothing, duplicate names, names changed, samples left out
    of a group, possible mixtures, groups without a table, HTML trees not written,
    and the groupings; the GUI reconciles what it staged against this instead of
    scraping HTML. `vsnp3_step2.py --version` (and step 1's, item 79) answers
    before the heavy imports: the GUI probes the version with a 10 s timeout and
    recorded an empty version when pandas took longer to import on a cold network
    file system.

## Packaging, wrappers and release checks (2026-10-16)

69. **`-r` is the reference everywhere the wrappers and the container README say
    so.** The four cluster wrappers (`internal/*_step1_cluster.sh`), the retired
    misc wrappers and `container/README.md` showed `-f *fasta -b *gbk` for the
    reference; `-f` is FASTA-to-FASTQ conversion, so a user following them
    converted the reference to reads and ran against a sourmash-chosen reference
    with no error. All now say `-r`, with a line explaining the difference. The
    step 1 help text is fixed with the step 1 items.
70. **The GPU wrapper no longer passes a flag the released code rejects.**
    `internal/git_gpu_vsnp3_step1.sh` added `--gpu`, which exists only on the
    parabricks-gpu-acceleration branch; on main argparse exited 2. It adds the
    flag only when `vsnp3_step1.py --help` lists it and otherwise says so and runs
    the CPU pipeline; it carries a `version="3.37"` line that `bump_version.sh`
    rewrites, prefers `vsnp3_parabricks_<version>` and falls back to the newest
    parabricks environment present. `docs/parabricks_gpu.md` says the flag is
    branch-only.
71. **The recipe sha256 cannot be forgotten again.** `internal/bump_version.sh`
    replaces the sha256 with the placeholder `REGENERATE_AFTER_TAG`, and
    `tests/test_release_consistency.py` fails until a real hash is written back
    (the 3.37 recipe sat with the 3.36 hash until the audit noticed). A 64-hex
    value is a WARN until `VSNP3_VERIFY_SHA256=1` downloads the tag tarball and
    compares it. `meta.yaml` carries the regenerate command above the hash.
72. **Pins.** `vcflib >=1.0.12` (was an exact pin that would break the solve when
    freebayes moves), `regex >=2023.0` (the old `2.5.110` was the module's internal
    version string, which conda never sees), `snp-dists` declared (the program
    `vsnp3_kernel_plots.py` runs; the script exits with a plain message when it is
    absent), `license: GPL-3.0-or-later` with `license_family: GPL3`. Both
    `meta.yaml` and `conda_build/vsnp3.yaml`; the release check forbids exact pins
    and non-SPDX licenses.
73. **Container: conda-forge base, working run script.** `container/Singularity.def`
    builds from `condaforge/miniforge3` with no `defaults` channel; the run script
    prints usage, runs step 1 when given options, or runs a named program (the old
    one called `vsnp3`, a command that does not exist). `container/README.md` links
    the definition instead of carrying a diverged copy, uses generic paths, and
    says `vsnp3:v3.37`; `bump_version.sh` and the release check cover it.
74. **Wrappers that could not run 3.37 are in `internal/attic/`** with a README
    saying why and what to use instead: the Azure step 1/2 wrappers and the misc
    single-sample worker and sourmash wrapper (Python 3.9 / pandas 1.3 / vsnp 3.14
    modules), the two misc callers of that worker (dead paths) and
    `phylo-snp_tree.sh` (its target script exists nowhere). `misc/` is empty and
    gone. The remaining cluster wrappers use a real here-doc comment instead of a
    triple-quote block that bash tried to execute, drop the undefined `$rflag`,
    and derive their environment hints from `${version}` (no more `vsnp3_3.26`,
    `vsnp3=3.31` or `vsnp3_vsnp3_3.37`).
75. **`bin/vsnp3_ksnp3_wrapper.sh`:** Inkscape 1.x export syntax with a fallback
    when inkscape is absent, a CPU count that never goes to zero, quoted file
    names, and a header saying it drives kSNP3 (the docs describe kSNP4).
76. **Documentation and hygiene.** README lists the programs the conda package does
    not install (picard, ncbi-datasets-cli, kraken2/Krona, the kSNP3 tool chain) as
    "install separately"; `.gitignore` typo fixed, a dead rule removed and the
    working-notes rule written as a glob; stray `.DS_Store` files and
    `bin/__pycache__` removed.
77. **Tests.** `tests/test_smoke.py` fails on a half-installed environment: nine
    recipe requirements were listed as optional, so a missing seaborn or plotly
    produced SKIP; only weasyprint (lazy, PDF only) and distro (not in the recipe)
    remain optional. `tests/test_release_consistency.py` checks
    `container/README.md`, fails on a stale wrapper `version=` line, warns on other
    version literals in wrappers, warns when `.github/workflows/ci.yml` is absent
    instead of skipping, warns on stray files, and builds its forbidden-word list
    at run time so the words appear neither in the source nor in its output.
    `tests/make_test_env.sh` offers a `recipe` leg (`conda_build/vsnp3.yaml` as
    declared) and falls back to it when ci.yml is absent. The local
    `vsnp3_3.27` environment lacks seaborn, so the smoke test is run in an
    environment holding every recipe requirement (118 of 118).

## Step 1 robustness (audit of 2026-10-16)

Checked on the SARS-CoV-2 and AF2122 test reads before and after, with and without
`VSNP3_MAX_CPUS=2`: VCF records, BAM alignments and every statistics field identical;
the step 1 HTML report differs only in the run-log echo of two renamed working files
(item 95). `tests/test_step1_units.py` (item 96) covers each fix.

78. **Step 1 opens every text file as UTF-8 with replacement.** Under `LC_ALL=C`,
    which every module sets, Python's default text encoding on Linux is ASCII, and
    freebayes and bcftools write the absolute reference path into the VCF header.
    An accent or en dash anywhere in the working directory or the reads path then
    raised `UnicodeDecodeError` after the alignment had completed: in the MQM
    rewrite and QUAL filter of `vsnp3_alignment_vcf.py`, the AC=1 count and header
    copy of `vsnp3_zero_coverage.py`, the nanopore writer of `vsnp3_longread.py`
    (which read the byte as U+FFFD and then could not write it), the validator's
    FASTQ and FASTA checks, the group reporter's VCF read and the run log itself.
    HPC jobs never saw it because `--export=NONE` unsets `LANG`; a workstation or
    the GUI inherits `en_US.UTF-8` and does. Every text-mode `open()` and
    `print(file=)` in the step 1 modules now says `encoding='utf-8',
    errors='replace'`, and `SeqIO.parse` is given an opened handle rather than a
    path, which also removes the two `ResourceWarning: unclosed file` lines. Step
    2's own reads were already fixed in 3.37 (item 67 for its remaining files).
79. **`--version` is answered before the heavy imports.** `vsnp3_step1.py`
    imported pandas, Biopython and dask before argparse saw `-v`, so on a slow
    network file system the vSNP GUI's provenance writer, which probes
    `vsnp3_step1.py --version` under a 10 s timeout, recorded an empty version.
    The flag is checked at the top of the file and the same text argparse printed
    (`vsnp3_step1.py: version 3.37`) is written in 0.01 s instead of 0.6 s.
80. **`-n/--SAMPLE_NAME` is honoured by the `_zc.vcf` name.** `Zero_Coverage`
    re-derived the sample name from the BAM file name, so `-n B26-0290_S3` named
    the BAM, the report and the statistics `B26-0290_S3` but the VCF `B26-0290`,
    and two forced names sharing a prefix collided in step 2; from the command line
    it also kept the directory part of the BAM path. The alignment passes the name
    it has; the regex on the BAM basename is the fallback, and the standalone
    `vsnp3_zero_coverage.py` takes `-n`.
81. **`-o output_dir` with relative paths runs.** The FASTQs were copied into the
    output directory and the process changed into it, after which a relative `-r
    ../ref/x.fasta` or `-b ../ref/x.gbk`, which validation had just accepted, was
    looked up from the new directory and `concat_fasta` raised `FileNotFoundError`.
    Every path argument is made absolute before the change of directory. A run
    without `-o` is unchanged.
82. **A missing or failing SPAdes no longer ends the run.** `-assemble_unmap` in an
    environment without spades.py raised `FileNotFoundError` inside the alignment,
    which `vsnp3_assembly.py` turned into `sys.exit(1)`: no VCF, report or
    statistics were written for the sample, and SPAdes' own output was discarded.
    The assembler's absence or failure is recorded ("failed assembly (SPAdes not
    found)" in the statistics, the reason in the run summary), the last lines
    SPAdes wrote go into the run log under a new "Tool output" heading, and the
    alignment carries on to the VCF as it already did for an assembly that
    produced no scaffolds.
83. **`krona_lca_all.py` parses the LCA mapping instead of destroying it.**
    `DataFrame.replace({' |:|': ''}, regex=True)` read the mate separator as the
    regex "a space, or a colon, or nothing" and deleted every space and colon from
    every mapping (`0:12 1763:5 |:| 0:10` became `01217635||010`), after which two
    `str.replace(..., regex=True)` calls raised `TypeError` on the first
    classified read with taxid 0, so the kraken2 wrapper exited at the Krona step.
    One plain replacement of the separator; the example resolves to taxid 1763.
84. **The help text no longer shows `-f` for the reference FASTA.** `-f` converts
    a FASTA into synthetic reads and analyses them; the reference is `-r`. The
    examples in `vsnp3_step1.py --help` said `-f *fasta -b *gbk`, so a user
    following them converted the reference genome to reads and aligned those to
    whichever reference sourmash chose, with no error. The examples say `-r`; the
    `-f` help text says what the flag is; `-f` combined with `-b` and no `-r` or
    `-t`, which is only ever this mistake, is refused by argparse with the
    explanation; the `-assemble_unmap` help in both files said "skip assembly" and
    now says the flag enables it; the conda line in the help names the current
    version instead of 3.13. (The cluster wrappers and the container README are
    item 69.)
85. **`vsnp3_fasta_to_fastq.py` gives each record the reads its own length calls
    for.** The whole-genome read count was generated for every record, so a
    two-chromosome Brucella came out at twice the requested coverage and a
    100-contig draft at a hundred times it. The random generator is seeded (0, or
    `--seed`), so the same FASTA gives the same reads every run; the read format
    is unchanged.
86. **seqkit statistics are read by column name.** `seqkit stats -a` was parsed by
    whitespace-split position, with Q20 and Q30 at indices 14 and 15, which
    assumes the `N50_num` column added in seqkit 2.4 while the recipe admits 2.1;
    on an older seqkit every sample read its mean quality as its Q30 and was
    Questionable. The default output is now parsed by header name, and a missing
    column is an error. The tab-separated form (`-T`) was tried first and
    rejected by the validation runs: seqkit 2.13 prints `Q20(%)` and `Q30(%)`
    there rounded to whole numbers while its default output keeps two decimals,
    so the workbook would have read "94%" where it had always read "93.61%".
    With the default output parsed by name every stored string is as before, on
    seqkit 2.10 and 2.13 alike.
87. **Spoligotype: at least one dask worker, and plain FASTQ accepted.**
    `int(cpus / 2)` was 0 on a machine with three cores or fewer, and
    `LocalCluster(n_workers=0)` is legal, so `client.gather()` waited for a worker
    that never came and the step hung for good without reaching its sequential
    fallback. Now `max(1, ...)`. `gzip.open` on an uncompressed `.fastq`, which
    the validator accepts, raised `BadGzipFile`; the opener is chosen by the gzip
    magic bytes.
88. **Re-running step 1 in the same folder works, and the file moves touch only
    this sample's files.** `shutil.move` raised when `alignment_<ref>/` already
    held a file of the same name, after the whole alignment had completed; the
    globs that gathered `*_nodup.bam`, `*_zc.vcf`, `*_annotated.vcf` and swept
    `*.sam`, `*.bai`, `*.fai`, `*_unfiltered_hapall.vcf` and the rest into `temp/`
    took every sample's files in the directory. Deliverables are moved by name,
    replacing what a previous run left; the cleanup patterns carry the sample name
    or the reference name. The final layout of `alignment_<ref>/` is unchanged.
89. **External tools run as argument lists.** The fifteen remaining `os.system`
    shell strings in `vsnp3_alignment_vcf.py` (faidx, bwa index and mem, minimap2,
    samtools fixmate/sort/markdup/index/fastq, freebayes-parallel, pigz) go through
    `vsnp3_run.run`, so `sample(1)_R1.fastq.gz`, a browser's duplicate download,
    aligns instead of giving `sh: syntax error` and a misleading "BWA SAM output is
    0 bytes", and a non-zero exit is reported where it happens. The command text
    recorded in the run summary is the text the shell string was; stderr goes
    where the old redirect sent it (the terminal, or `/dev/null`) unless it is
    kept. freebayes-parallel hands its arguments to GNU parallel, which runs them
    through a shell again, so the paths it receives are shell-quoted, a no-op for
    ordinary names.
90. **Thread counts follow `VSNP3_MAX_CPUS`.** Every aligner, samtools, bcftools
    and freebayes-parallel call said 8 threads whatever the machine or allocation.
    The alignment uses one count, `min(8, worker_count())`: 8 on a large machine
    as always, fewer under the cap the GUI sets or on a small machine. bwa mem's
    output depends on its thread count (its batch is 10 Mbp per thread and the
    insert-size distribution is estimated per batch; 160 of 600,000 read pairs
    moved between `-t 8` and `-t 2`), so when the count is not 8 the batch is
    pinned to the 8-thread size with `-K 80000000`, which was verified to give the
    same SAM, VCF and statistics as an uncapped run. The kraken2 wrapper's
    `cpu_count() - 2` is floored at 1.
91. **freebayes-parallel region failures are visible.** Its exit status is
    vcfuniq's, so a region whose freebayes was killed left a VCF with records
    elsewhere and nothing there, which step 2 reads as "matches the reference".
    Its stderr is kept in the run log, and after the call every contig with mapped
    reads (samtools idxstats) is checked for at least one record; a contig without
    one is named in a NOTE in the summary and the report, not failed on.
92. **Step 1's group lookup reports, and never ends the run.** The lookup caught
    four exception types; `BadZipFile`, `ImportError`, `AttributeError` or
    `AmbiguousReferenceError` from the reporter killed the sample after its
    `_zc.vcf` existed. And a reference type with no defining-SNP workbook was
    `pd.ExcelFile(None)`, reported as "group lookup FAILED: Invalid file path or
    buffer object type". The reporter raises `NoDefiningSnpWorkbook`, which step 1
    reports as "group file not provided"; every other exception from this
    reporting-only step is caught and reported as a failed lookup. The step 1
    group rule is step 2's: MQ at or above the threshold (it was strictly above)
    and a single-base ALT (it was unchecked), both from `vsnp3_version`.
93. **`vsnp3_run.run_pipe` drains the first process's stderr on a thread.** It was
    read only after the second process had finished, so once bcftools mpileup had
    written about 64 KB of warnings (many contigs, read-group notes) it blocked on
    stderr, stopped feeding stdout, and the two processes waited on each other for
    good. A consumer that fails to start no longer leaves the producer writing
    into a pipe nobody reads.
94. **Brucella MLST: picard checked up front, unreadable VCF reported, cleanup
    limited to the run's own files.** picard is not in the vsnp3 recipe; without
    it the pipeline ran four more commands and died on `df['POS']` with a bare
    `KeyError`. The required programs are checked before anything is written and
    named in the message; `read_vcf` raises with the file name and cause instead of
    returning a column-less frame (a header-only VCF is still a valid empty
    result); the `ST1-MLST*`, `*-mlst*` and `*_mlst.vcf.idx` globs, which deleted
    anything in the directory with those substrings, are replaced by the explicit
    list of files the run writes.
95. **Small fixes.** `vsnp3_zero_coverage.py` skipped header lines with
    `read_csv(comment='#')`, which also cut any record short at a `#` in a sample
    or contig name; the header is skipped by line and reused for the `_zc.vcf` (one
    read instead of two). The concatenated reference is written under
    `<stems>.vsnp3_reference.fasta` when `<stems>.fasta` already exists in the
    working directory, so a user's own `-r` FASTA is neither rewritten nor moved
    into `alignment_<ref>/` (the second sample run in that folder used to find no
    reference); the stem before the first `.` is the same, so `alignment_<ref>/`
    and the statistics do not change. The QUAL gate on freebayes output imports
    `LONGREAD_CALLER_MIN_QUAL` instead of a literal 20. The working files
    `v_header.csv`, `v_annotated_body.csv`, `markduplicate_stats.txt`,
    `chrom_ranges.txt` and `temp_fastq_seqkit_stats.txt` carry the sample name, so
    two samples run by hand in one folder cannot read each other's header,
    duplicate statistics or regions (the two file names appear in the run log's
    command lines, which is the one visible change in an ordinary run).
    `vsnp3_download_GCA_fasta_get_metadata.py` reports exhausted retries and exits
    1 instead of an `AttributeError` on `acc`; `vsnp3_alignment_vcf.py` requires
    `-r`; `vcf_search.sh` splits with `tr` and `awk` rather than GNU-only `\n` and
    `\t` in sed.
96. **Tests.** New `tests/test_step1_units.py`, 86 checks in 15 sections, pure
    Python with external programs mocked: the krona mapping, seqkit tables of two
    layouts, the sample name and VCF read of zero_coverage, dask worker count and
    FASTQ opener, `run_pipe` with 300 KB of stderr, SPAdes absent and failing,
    per-record read counts and seeding, the group rule at the MQ boundary and with
    a multi-base ALT, thread cap and bwa batch pin, replacing moves, the reference
    name, contigs without records, the step 1 command line (`--version` before the
    heavy imports, help text, `-f` with `-b` refused), bruc_mlst, kraken2 threads,
    exhausted downloads, non-ASCII FASTA/FASTQ, `vcf_search.sh`.

## Step 2 summary trimmed, mixed-sample flag stricter (2026-10-16)

97. **Three summary sections removed.** The step 2 summary no longer shows the
    "Remove list entries that matched no VCF file" box, the "Possible mixed
    samples" section or the "Calls read from inside multi-base records" section.
    The remove list is usually the reference type's and names samples of other
    runs too, so its box was mostly noise. All three are still printed in the
    terminal at the end of the run, and the accounting JSON still carries
    `remove_list_entries_unmatched` and `possible_mixed_samples`, so the vSNP GUI
    sees no change.
98. **Possible mixed samples are marked in the Groupings table.** A flagged
    sample's name is highlighted (bold, red on cream), and after its groups the
    row adds "N confident mixed calls against M passing SNPs". One line above the
    table explains the highlight, shown only when a sample is flagged. The step 1
    notice and the terminal message point to the Groupings table.
99. **Only the clearest mixtures are flagged.** `vsnp3_version.mixed_signature`
    now needs at least 20 confident mixed calls (`MIXED_MIN_AC1`, was 10) and
    more than 5 per passing SNP (`MIXED_AC1_RATIO`, was 1.0). On a B. suis run the
    old rule flagged 13 samples; 7 had 1.0 to 2.8 calls per passing SNP and are no
    longer flagged, and the 6 at 5.5 to 13.6 still are. Validation limits: clonal
    samples had at most 0.41 calls per SNP and known mixtures at least 20, so 5 is
    12 times the first and a quarter of the second. The same rule drives the
    step 1 notice. The rule for samples with no passing SNP (every call mixed)
    is unchanged.
100. **The GUI marker is back in `vsnp3_step2.py`.** Item 51 says both pool files
    name `VSNP3_MAX_CPUS`, which is what the vSNP GUI's `apply.sh` checks before it
    rewrites the pool size. `vsnp3_step2.py` had lost the name when its pool
    moved to `worker_count()`. It is now in a comment at that call again.
101. **Tests.** `tests/test_step2_inputs.py` grows from 186 to 191 checks: the
    rule's new limits, the 13 B. suis samples (6 kept, 7 dropped), and a new
    section [18] that builds the summary and checks that the three sections are
    gone and that a mixed sample is highlighted with its counts after its groups.

## Files changed

- `bin/vsnp3_step2.py` — reader, accounting, summary (rewritten); eight-column read and blank-line count (items 28, 30), tree mode and HTML-tree failures in the summary (32, 33); CIGAR kept for multi-base records, multi-base and mixed-sample sections in the summary (38, 39); platform read from the header, long-read samples in the summary (44); "Not analysed" list moved to the end, `removed_samples` gone (49); `worker_count` (50); `single_base` column, `DefiningSnpsError` exit, deletion at the end, summary boxes and `group_failures`, accounting JSON, `--version` before imports, thresholds from `vsnp3_version` (52-68)
- `bin/vsnp3_input_validator.py` — validation (rewritten); strict majority (29); truncated files (40)
- `bin/vsnp3_group_on_defining_snps.py` — metadata, no-SNP samples, phases, timeout wrapper; largest group first, frame handed to html_tree; `html_tree_build` and failure list (33); `_embedded_calls`, `_allele_bases` (38), `_confident_mixed`, `mixed_samples` (39); `worker_count`, GUI marker comment (50, 51); `read_metadata` blank rows, `safe_fallback_name`, `_check_defining_snps`, `DefiningSnpsError`, one `_passing` pass with `_single_base`, mask-based parsimony test, `undetermined_in_group`, `rewrite_errors`, `_html_tree_timeout` (54-61)
- `bin/vsnp3_html_tree.py` — `SnpMatrix` tree builder, bitmask clade keys, batched plotly calls, `df=` parameter, recursion limit; `read_table`, `raise_errors` (33, 34); `NON_SAMPLE_ROWS` (53), iterative `count_leaves` (61), one plotly bundle (65)
- `bin/vsnp3_fasta_to_snps_table.py` — RAxML slots and version cache, `snp_column_stats`, direct workbook writer, `as_read_back`; column A format (27), names column (34), `bootstrap_replicates` and comments (31); unused `Tree.cpu_count` removed (50); `leaf_order` via Bio.Phylo and the alignment-vs-tree check (52), `html_tree_skipped` (59), formula-safe names (66)
- `tests/regression/bench_tree_phase.py` — new; `tests/regression/README.md` describes it; compares with `read_table` (37)
- `bin/vsnp3_file_setup.py` — `PhaseTimer`, `format_elapsed`, usability verdicts, `HTML_ONLY_KEYS`, report tile note; SNP caller and long-read rows in the report (46); `worker_count` (50)
- `bin/vsnp3_step1.py` — verdicts with reasons, no-SNP notice (35), mixed-sample notice (39); `-p` help (46)
- `bin/vsnp3_zero_coverage.py` — `step2_snp_counts`, `no_snp_notice` (35), `confident_mixed_count`, `mixed_sample_notice` (39)
- `bin/vsnp3_alignment_vcf.py`, `bin/vsnp3_run.py` — header-only freebayes VCF accepted, `require_output(vcf_header=True)` (36); nanopore path: minimap2 options and read group, no duplicate marking, flagstat counts, `vsnp3_longread` calling, unmapped reads, `SNP Caller` and long-read statistics (42-46)
- `bin/vsnp3_longread.py` — new: the nanopore caller and the evidence QUAL, with a command line (42, 43, 47)
- `bin/vsnp3_version.py` — usability limits; version 3.37 (26); `mixed_signature` and its limits (39); `LONGREAD_*` limits, `NANOPORE_QUAL_OFFSET` removed (42, 43)
- `bin/vsnp3_reference_options.py` — basename matching; exact header-ID detection and `detection_note`, UTF-8 reads (63, 67)
- `bin/vsnp3_annotation.py`, `bin/vsnp3_path_adder.py` — UTF-8 handles (67)
- `bin/vsnp3_excel_merge_defining_snps.py` — basename, cell checks (64)
- `bin/vsnp3_vcf_merge_to_fasta.py` — thresholds imported, `_zc` suffix, single-contig check (66)
- `bin/vsnp3_vcf_annotation.py` — temp file names, output name (66)
- `bin/vsnp3_remove_from_analysis.py` — `entries` kept for reporting (62)
- `README.md`, `docs/instructions/additional_tools.md` — viewer link
- `conda_build/meta.yaml`, `container/Singularity.def`, `container/README.md`, `internal/*.sh` — version 3.37 (26)
- `LICENSE`, `.gitignore`, `internal/bump_version.sh` — added (26)
- `tests/test_step2_inputs.py` — new; D-02 to D-12 checks (37); L-2, L-3, D-06 checks (41); long-read section [16] (48)
- `README.md`, `docs/instructions/additional_tools.md` — nanopore reads (42-47)
- `tests/test_release_consistency.py` — works outside a git checkout (26)
- `conda_build/meta.yaml`, `conda_build/vsnp3.yaml` — vcflib floor, regex floor, snp-dists, SPDX license, sha256 comment (71, 72)
- `conda_build/instructions` — placeholder flow, container/README.md, post-release GUI retest note (71)
- `internal/bump_version.sh` — sha256 placeholder, container/README.md (71, 73)
- `internal/git_gpu_vsnp3_step1.sh` — `--gpu` guard, `version=` line, environment fallback (70)
- `internal/vsnp3_step1_cluster.sh`, `git_`/`nbaf_`/`scomp_vsnp3_step1_cluster.sh` — here-doc comment, `-r`, `$rflag` removed (69, 74)
- `internal/vsnp3_step1.sh`, `vsnp3_step2.sh`, `scomp_vsnp3_step1.sh`, `scomp_vsnp3_step2.sh`, `scomp_git_vsnp3_step2.sh`, `git_vsnp3_step1.sh`, `git_vsnp3_step2.sh` — hints derived from `${version}` (74)
- `internal/attic/README.md` — new; `internal/attic/*.sh` — moved from `internal/` and `misc/` (74)
- `container/Singularity.def`, `container/README.md` — miniforge base, run script, link, generic paths, 3.37 tags (69, 73)
- `bin/vsnp3_ksnp3_wrapper.sh` — Inkscape 1.x, CPU guard, quoting, kSNP4 note (75)
- `bin/vsnp3_kernel_plots.py` — snp-dists check (72)
- `docs/instructions/additional_tools.md` — kSNP3 note (76); `docs/parabricks_gpu.md` — branch-only flag (70)
- `.gitignore` — typo, dead rule, `*-handoff.txt` glob (76, 77)
- `tests/test_smoke.py`, `tests/test_release_consistency.py`, `tests/make_test_env.sh` — (77)
- `tests/test_step2_inputs.py` — section [17], robustness fixes (52-66)
- `bin/vsnp3_step1.py` — `--version` before the imports (79); UTF-8 opens and run log (78); absolute paths with `-o` (81); help text, `-f` with `-b` refused (84); group lookup outcomes and catch-all (92); "Tool output" in the run log (82, 91); reference name helper (95)
- `bin/vsnp3_alignment_vcf.py` — argv command execution with identical recorded text (89); `step1_threads`, bwa `-K` pin (90); freebayes-parallel stderr and contig check (91); `move_replacing`, named deliverables, narrowed cleanup (88); SPAdes failure path (82); `concatenated_reference_name`, `LONGREAD_CALLER_MIN_QUAL`, sample-named working files (95); `-r` required, `-assemble_unmap` help (84, 95); UTF-8 opens (78)
- `bin/vsnp3_zero_coverage.py` — `sample_name_for`, `sample_name=` and `-n` (80); `read_vcf_lines`/`vcf_body_frame` (95); sample-named temp files (95); UTF-8 opens (78)
- `bin/vsnp3_assembly.py` — `spades_argv`, `_run_spades`: no `sys.exit`, log tail kept (82)
- `bin/vsnp3_run.py` — `stderr=` on `run` (89); threaded drain in `run_pipe` (93)
- `bin/vsnp3_fastq_stats_seqkit.py` — `seqkit stats -a -T` by header name, `seqkit_display` (86); sample-named temp file (95)
- `bin/vsnp3_group_reporter.py` — `NoDefiningSnpWorkbook`; step 2 record rule (92); UTF-8 read (78)
- `bin/vsnp3_spoligotype.py` — `half_cpus`, `open_fastq_text` (87)
- `bin/vsnp3_bruc_mlst.py` — tool check, `read_vcf` raises, explicit cleanup (94)
- `bin/vsnp3_kraken2_wrapper.py` — threads floored at 1 (90); UTF-8 opens (78)
- `bin/krona_lca_all.py` — `best_taxid`, regex replace removed (83)
- `bin/vsnp3_fasta_to_fastq.py` — `reads_for_length` per record, `--seed` (85)
- `bin/vsnp3_longread.py`, `bin/vsnp3_input_validator.py` — UTF-8 opens (78)
- `bin/vsnp3_download_GCA_fasta_get_metadata.py` — exhausted retries reported, exit 1 (95)
- `bin/vcf_search.sh` — portable splitting (95)
- `tests/test_step1_units.py` — new (96)
- `bin/vsnp3_step2.py` — remove-list box, mixed-sample and multi-base sections out of the summary; mixed samples highlighted in the Groupings table; terminal wording; `VSNP3_MAX_CPUS` comment (97, 98, 100)
- `bin/vsnp3_version.py` — `MIXED_MIN_AC1` 20, `MIXED_AC1_RATIO` 5.0 (99)
- `bin/vsnp3_zero_coverage.py` — step 1 mixed notice wording (98, 99)
- `bin/vsnp3_step1.py`, `bin/vsnp3_group_on_defining_snps.py` — comments on the mixed rule (99)
- `tests/test_step2_inputs.py` — mixed-rule limits, section [18] (101)
- `CHANGELOG_3.37.md` — this file
