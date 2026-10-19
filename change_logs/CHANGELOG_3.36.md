# vSNP3 3.36 — release tag summary

Concise list of what changed. Detail follows below.

**Corrected results — re-analysis needed**
1. Minus-strand amino acid calls were wrong in 3.16/3.34/3.35 (0% correct on the minus strand; ~50% of bacterial CDS). Affects step 2 SNP tables **and** step 1 `*_annotated.vcf`. The `_zc.vcf` calls themselves were never wrong, so no re-alignment is needed.
2. HTML trees reported "0 SNPs" for 9 of 23 shipped references (any non-`NC_` contig: orygis, mtbc0, LT708304, para-CP033688, 4 Brucella).
3. A chromosome missing from the gbk ended step 2 mid-run with exit status 0.
4. Illumina samples with one long read were routed to the nanopore caller and had +100 added to every QUAL.
5. Missing QUAL was treated as 0, which rewrote the position to the **reference base** instead of a no-call.
6. `AC=1` sites with an unrecognized REF/ALT were published as homozygous, and the first one aborted all remaining het sites in that sample.
7. A tree sample missing from the alignment was replaced by a substring match — one isolate's SNPs under another's name.

**Data no longer destroyed**
8. `vsnp3_step2.py -wd <dir>` deleted the VCFs in `<dir>` even when `<dir>` was a shared database, and `-k` never worked. Files this run did not create are now never deleted; `-k` works.
9. Input VCFs were rewritten in place, corrupting a trailing INFO value to 999 and forcing genotypes to `1/1`.
10. `vsnp3_group_reporter.py` deleted the sample's `_zc.vcf` on a parse error and reported it as "group file not provided".
11. A failed zero-coverage merge silently emitted the unmodified VCF as `_zc.vcf`.

**Reproducibility**
12. Step 2 output is now byte-identical across runs: seeded bootstrap, deterministic multi-allele tie-break, sorted globs.

**Fixed tools**
13. `vsnp3_excel_merge_defining_snps.py` could never complete a run (malformed format string).
14. `vsnp3_group_on_defining_snps.py -t` raised `AttributeError` immediately; standalone also applied no filters.
15. `vsnp3_kraken2_wrapper.py --help/--version` exited 1.
16. Every default step 1 run hit an `UnboundLocalError` while reporting program versions; step 1 hid it by matching the message text.
17. A step 2 group in which every position carries the same call in every sample killed the whole run with `KeyError: 'abs_pos'`, naming one sample as if its VCF were malformed. Such a group has nothing to align; it is now reported and skipped, and every other group still produces its tree.
18. A VCF matching several rows of the metadata worksheet is now reported. The lookup takes the first matching row, so the others were discarded silently -- and because the file still ends up with one name, no count changed and nothing anywhere mentioned it. Reported only when the rows disagree.
19. Two VCFs resolving to one sample name no longer aborts the run. Names are now resolved in their own pass before any sample is processed, one file per name is analysed -- the last of the colliding files in sorted order, chosen explicitly rather than left to directory order -- and the HTML summary carries a note naming the file used and the ones skipped, so the duplicate can be corrected in the metadata worksheet afterwards. The old guard was fatal and its message listed every input file rather than the colliding ones: it built a `Counter` of the basenames and printed its keys, never its counts, which on a 7,700-sample run was a 167 KB message identifying nothing.
20. Records with no MQ or QUAL are no longer reported. A zero-coverage record carries neither by construction, so in a `_zc` VCF this is most of the file and says nothing about the sample; at 7,992 samples it was 7,700 lines of log burying the traceback beneath it. They are still excluded from position selection rather than read as zero, and still reach the table as `-`.

**New**
21. Step 1 writes the report twice: a self-contained HTML file, and a PDF rendered from that same HTML by WeasyPrint.  LaTeX is gone, and with it the texlive, svgwrite and cairosvg dependencies; `weasyprint` and `jinja2` are declared in their place.  Nothing the LaTeX report showed was lost — the sourmash similarity ranking, its one part that never reached the HTML, is now a report section.
22. New flags: `--assume_gt_only_quality`, `--first_sample_only`, `--bootstrap_seed`, `-n/--no_filters`, `vsnp3_version.py --check`. **No option was removed** — every 3.35 command line still works.
23. Reports disclose no absolute paths.  The step 1 run log recorded each command line verbatim, including the full path of the reference and the reads, and the step 2 summary printed the `--gbk`, `--metadata`, `--defining_snps` and `--remove_by_name` paths as given -- so a report sent outside the lab described the sending machine's directory layout and account name.  Filenames are kept; directories are not.
24. Whitespace in any path is refused up front instead of silently producing 0× coverage.
25. `NCBI_EMAIL` now required by the accession downloader (was a hardcoded fake address).

**Performance — no result changes; each verified data-identical**
21. Step 1 ~60 s faster per sample: one `conda list --json` instead of 22 calls (9.85 s → 0.60 s), and a `FASTQ_Stats` run whose result was discarded is gone (4 redundant full FASTQ passes).
22. Coverage stage memory 1,482 MB → 17 MB per sample: `samtools depth -a` streamed into numpy arrays instead of four whole-genome dicts and an outer merge on a 4.35 M string index. Matters most on WSL and laptops.
23. Density filter O(P²) → linear sweep: 15.94 s → 0.007 s at 40,000 positions.
24. `--find_new_filters` 112× faster and no longer runs twice per group.
25. HTML tree builder 2.4× faster (per-round SNP-set memoization), tie-break preserved. Still O(S³) — see below.
26. RAxML binary resolved once via `shutil.which` instead of executing up to 7 candidates per group.

**Packaging / infrastructure**
27. `meta.yaml` 3.23 → 3.36; declared `pyyaml`, `jinja2`, `matplotlib-base`, `seaborn`. **Rebuild the env, don't just update the package.**
28. Removed a committed 1.4 GB container image — source tarball is now 4.9 MB.
29. One version source (`bin/vsnp3_version.py`) replacing 29 literals; `internal/bump_version.sh`; `vsnp3_version.py --check` verifies a whole install.
30. First test suite and CI: 94 smoke + 39 annotation checks, release-consistency gate, a run-twice determinism gate in `tests/regression/`, and a conda build of the recipe as written.
31. `bin/get_stats.sh` and `bin/organize_directories.sh` could `rm -r` in the wrong directory after a failed `cd`; fixed before `build.sh` began installing them.

**Known limits**
32. Tree builder is still O(S³): ~235 s per group at 200 samples against a 900 s timeout. Pairwise-count cache deferred to 3.37.
33. External tool thread counts are still hardcoded (up to 3.3× oversubscribed); `Setup.cpus` is computed and never used. Deferred to 3.37.
34. Sample names still truncate at the first `_` or `.`; collisions are now detected and refused, but the truncation itself needs the shell wrappers changed in step. Deferred to 3.37.
35. **Outstanding before tagging:** add a `LICENSE` file (declared in `meta.yaml`, has never existed) and regenerate the `sha256`.

---

# vSNP3 3.36 — detail

A correctness release. It began as one reported defect — minus-strand amino acid
calls — and a review of the rest of the pipeline found further defects that
silently produced wrong data, destroyed input files, or made output depend on
directory order. **Results will change relative to 3.35, and some prior results
should be re-analysed.** Sections 1 and 2 say which.

---

## 1. Re-analyse these

### Every minus-strand amino acid call from 3.16, 3.34 and 3.35 is wrong

Stock 3.35 reported the plus-strand reading of every minus-strand codon. Measured
against each reference's own `/translation` qualifier on H37Rv single-CDS
positions:

| | plus strand (n=336) | minus strand (n=341) |
|---|---|---|
| 3.35 | 100.0% correct | **0.0% correct** |
| 3.36 | 100.0% | 100.0% |

Half of H37Rv's 3906 CDS are on the minus strand. The canonical control,
katG Rv1908c `NC_000962.3:2155168 C>G`, was reported as **A315G** and is
**S315T**. On bacterial references roughly 30% of in-feature SNPs change
synonymous↔nonsynonymous.

**Both step 1 and step 2 outputs are affected.** Verified by running the 3.35 and
3.36 annotators over the same VCF against the AF2122 reference — 517 real
single-base coding SNPs:

| | minus strand (n=175) | plus strand (n=136) |
|---|---|---|
| `snp_aa` changed | **100%** | 1.5% |
| synonymous↔nonsynonymous flipped | **58.3%** | 0.7% |

Example, gene `ino1`: 3.35's annotated VCF said `CTT→TTT L213F nonsynonymous`;
the correct call is `AAG→AAA K213K silent`. A different biological conclusion.

So the step 1 `*_annotated.vcf` files carry the same wrong amino acid calls as the
step 2 tables — they contain `ref_aa`, `snp_aa`, `mutation_type` and both codons,
all of which the defect corrupted.

**What does and does not need redoing:**

- **No re-alignment or re-calling.** The `_zc.vcf` files — positions, REF, ALT,
  QUAL, MQ, DP — were never affected. Keep them.
- **Re-generate the annotated VCFs** from the existing `_zc.vcf` with
  `vsnp3_vcf_annotation.py -b <ref>.gbk -c <sample>_zc.vcf`.
- **Re-run step 2** on the existing `_zc.vcf` files to get corrected SNP tables.

Pinning 3.34 is **not** a workaround — the regression entered during the 3.34 fix
series, so 3.34 is affected too.

To check any install, run the katG control rather than trusting the version string:

```bash
vsnp3_annotation.py -g /path/to/Mycobacterium_H37/NC_000962.gbk \
    -p NC_000962.3:2155168 -s G
```

A correct install reports `AGC` → `ACC`, `S315T`, `p.Ser315Thr`. A defective one
reports `GCT` → `GGT`, `A315G`.

(An earlier draft of these notes suggested `grep -c reverse_complement
vsnp3_annotation.py`. That is wrong and would mislead: 3.36 takes the codon from
`feature.extract()`, which is already strand-aware, so it needs no explicit
reverse-complement call and greps to 0 exactly like broken 3.35 does. Only the
functional check distinguishes them.)

### Trees from nine of the shipped references reported zero SNPs

The HTML tree kept only table columns whose name began with `NC_`, so any
reference whose contigs are not RefSeq `NC_` accessions produced a star tree
reporting "Total SNPs in tree: 0" — written out as a normal result. Affected:

| reference | contigs |
|---|---|
| `Mycobacterium_orygis` | CP063804 |
| `mtbc0_v1.1` | MTBC0 |
| `LT708304` | LT708304 |
| `para-CP033688` | CP033688 |
| `Brucella_abortus3` | NZ_CP007682, NZ_CP007683 |
| `Brucella_melitensis-bv1b` | NZ_CP018508, NZ_CP018509 |
| `Brucella_melitensis-bv3` | NZ_CP007760.1, NZ_CP007761.1 |
| `Brucella_suis3` | NZ_CP007719, NZ_CP007718 |
| `Brucella_neotomae` | KN046827 |

AF2122 and H37Rv are both `NC_`, which is why this went unnoticed. Any HTML tree
published from one of the references above is meaningless and should be re-run.

### Runs that ended early with a success status

A chromosome with no GenBank record made step 2 call `sys.exit(0)` mid-run: no
tables, no trees, exit status 0, and nothing for a wrapper or scheduler to
detect. If you have runs whose output directory is missing its tables, this is
the likely cause. Those positions are now reported as unannotated and the run
continues.

### Illumina samples routed through the nanopore caller

Platform detection tested the **maximum** read length, so a single long read in
an Illumina library sent the whole sample through minimap2 and bcftools *and*
added +100 to every QUAL — making those values incomparable to the thresholds
applied downstream. Detection now uses average read length, and the presence of
an R2 file settles it. Affected samples must be re-run. The decision and the
statistic behind it are now written to the run log and the Excel stats.

### Input VCFs were modified and deleted

- `vsnp3_step2.py -wd <dir>` **deleted the VCF files in `<dir>`**. On a shared
  reference database that removed everyone's inputs; the only copy left was
  `vcf_starting_files.zip` in the output directory. `-k/--keep_ind_vcfs` existed
  to prevent this and was parsed but never used.

  3.36 keeps the useful half of that behaviour and drops the dangerous half. The
  files are still archived to `vcf_starting_files.zip` and the loose copies are
  still cleaned up by default, because leaving them clutters the working
  directory. But only files this run owns are ever deleted — either copies it
  made itself under `--output_dir`, or files inside the working directory it was
  pointed at. A `--wd` pointing somewhere else is left alone with an explanation.
  `-k` suppresses deletion entirely.

  | invocation | files live in | deleted? |
  |---|---|---|
  | `-wd .` | the working directory | yes |
  | `-wd /db -o out` | copies in `out/` | the copies; `/db` untouched |
  | `-wd /db` | `/db` | no — reported, not deleted |
- Normalization rewrote each input file in place. For records whose `FORMAT` was
  exactly `GT`, it overwrote the **last numeric INFO value** with 999 — usually
  `MQ`, which then passed the mapping-quality filter unconditionally — and forced
  every genotype to `1/1`. A `0/1` heterozygous record became `AC=999` and `1/1`,
  which then failed both the `AC==2` and `AC==1` tests and vanished from the
  analysis.
- Multi-sample VCFs were silently reduced to their first sample.

If you have VCFs that passed through 3.35 or earlier, prefer the pristine copies
in each run's `vcf_starting_files.zip`.

### `_zc.vcf` files that silently lost their zero-coverage information

If inserting zero-coverage positions failed, step 1 fell back to copying the
unmodified VCF as `_zc.vcf`. That file is precisely what distinguishes "no
coverage here" from "matches the reference here", so the fallback asserted a
reference call at every uncovered position with no marker. It now raises and
writes nothing.

---

## 2. Other output changes

| change | effect |
|---|---|
| Sample and position ordering is now deterministic | Row and column order in tables and FASTA records follows VCF filename order rather than directory order. Calls are unchanged. |
| Amino acid calls at multi-allele positions | Previously the ALT that reached the annotator was whichever VCF was read last, so identical input could give different published calls. Now the most frequent ALT, ties broken alphabetically, with a `multi-allelic [G,T]` marker. |
| Bootstrap support values | Were unseeded: two runs of the same data gave one clade 100% and then 98%. Now seeded, so they are reproducible — and will differ from any prior run. |
| `AC=1` positions with an unrecognized REF/ALT pair | Were published as confident homozygous calls, and the first such row abandoned every remaining heterozygous position in that sample. Now `N`, and reported. |
| Missing QUAL | Was filled with 0, which put the position below the N threshold and rewrote it to the **reference base** — turning an absent measurement into positive evidence of a reference match. Now a no-call (`N`). |
| Missing MQ | Was filled with 0 and silently failed the filter. Now reported and excluded; a VCF with no MQ field at all excludes the sample with an explanation instead of contributing nothing. |
| A tree sample absent from the alignment | Was replaced by a substring match, so `19-01` matched `19-011` and one isolate's SNP row was written under another isolate's name. Now an error. |
| Multi-chunk tables with `--show_groups` | Conditional formatting was shifted one column on chunks 2..N. |
| Deterministic GenBank selection | When two gbk files define the same chromosome, the winner is now stable and a warning names both files. |
| Step 1 report format | Was a PDF via `pdflatex`, which was invoked unchecked and with texlive undeclared, so a machine without LaTeX silently produced no report -- and the R2 average read length it printed was R1's value. Now an HTML file rendered from the same stats dict that produces the workbook, plus a PDF rendered from that HTML, so the two cannot disagree. Both results are checked, and a missing or failing PDF renderer is reported rather than swallowed. |
| Error exits are non-zero | Thirteen error paths exited 0. Wrappers that treat exit 0 as success may start reporting failures that were previously invisible — that is the intent. |
| Sample-name collisions are refused | Names truncate at the first `_` or `.`, so `TB_2024_001` and `TB_2024_002` both became `TB` and one sample silently vanished. Now detected and refused. Changing the truncation itself is deferred to 3.37, because the shell wrappers do the same truncation and users' metadata files are keyed on the truncated name. |

---

## 3. Command line

| flag | change |
|---|---|
| `-k`, `--keep_ind_vcfs` | Now functional. It was parsed and never read in 3.35. Pass it to keep the VCF files that were read; by default they are deleted once archived, since leaving them clutters the working directory. |
| `-f`, `--fix_vcfs` | Now writes normalized **copies** to `vcf_normalized/` and leaves the inputs alone. Reports files it could not normalize and exits non-zero. |
| `--assume_gt_only_quality` | New. For assembly- or consensus-derived VCFs whose `FORMAT` is `GT` only and which therefore carry no quality data: substitutes QUAL, leaves genotypes and INFO untouched, and records `##vsnp3_synthesized_quality` in the header. Replaces the previous silent, and incorrect, substitution. |
| `--first_sample_only` | New. Required to analyse a multi-sample VCF, which is otherwise refused rather than silently reduced to sample one. |
| `--bootstrap_seed` | New on `vsnp3_html_tree.py`. |
| `-n`, `--no_filters` | New on `vsnp3_group_on_defining_snps.py`, defaulting to False to match `vsnp3_step2.py`. Run standalone, it previously applied **no** position filters, so the same inputs gave a different matrix depending on the entry point. |
| Whitespace in paths | Refused up front. The aligner and caller command lines interpolate paths into shell strings unquoted, so a space produced an empty SAM, an empty BAM, an empty VCF, and a run that reported success at 0× coverage. |
| `NCBI_EMAIL` | `vsnp3_download_fasta_gbk_gff_by_acc.py` now requires this environment variable. It previously sent `mickey_mouse@gmail.com`, which violates NCBI's E-utilities policy and risks throttling every user on a shared institutional address. `NCBI_API_KEY` is honoured if set. |

Two entry points that could not complete a single run are fixed:
`vsnp3_excel_merge_defining_snps.py` built both of its output filenames with a
malformed format string and died at write time after doing all its work; and
`vsnp3_group_on_defining_snps.py -t <reference>` raised `AttributeError`
immediately.

---

## 4. Testing

There was no CI, and the only automated test in the release path was
`vsnp3_step1.py -h`. That is why a defect making every minus-strand call wrong
shipped in two consecutive releases.

- `tests/test_annotation.py` — 39 checks: synthetic plus/minus/spliced/ragged/
  tRNA/pseudogene/overlapping/non-ATG cases, the canonical drug-resistance
  controls, and a `/translation` oracle across five shipped references. Moved out
  of `bin/`, where `build.sh`'s `mv bin/*.py` had been shipping it to every user.
  Set `VSNP3_REQUIRE_DATA=1` to make a missing dataset fail rather than skip.
- `tests/test_smoke.py` — 93 checks needing no external data: every module
  compiles with no `SyntaxWarning`, every literal format string parses, every
  script's `--help`/`--version` exits 0 and all 26 agree on the version.
- `tests/test_release_consistency.py` — version agreement across the recipe, the
  READMEs, the container definition and the wrappers; declared dependencies; and
  nothing oversized in the source tarball.
- `.github/workflows/ci.yml` — unit and release-consistency on every push;
  integration, including a run-twice byte-comparison, on pull requests.
- `Annotation.self_test()` translates the codon the module actually reports and
  compares it against `/translation`. It is proven to fire: reintroducing the
  plus-strand slice produces
  `### ANNOTATION SELF TEST FAILED (28 of 60 sampled codons)`.

**What this does and does not cover.** Be clear-eyed about it: the suite covers the
annotation engine thoroughly, and covers every script's import graph, argparse
wiring and version reporting. It does **not** functionally test most options — of
86 distinct command-line options across the suite, only `--help`/`--version` are
exercised everywhere, plus the one step 2 path (`-wd`, `-t`, `-html_tree`, `-a`)
and step 1 (`-r1`, `-r2`, `-t`) used by the end-to-end tests. Options with no
functional test include `--filter_density`, `--hash_groups`, `--show_groups`,
`--remove_by_name`, `--all_vcf`, `--assemble_unmap`, `--spoligo` and the
`-w/-x/-y` thresholds. Nor is this pytest: the suites are plain scripts with a
`check()` helper. Converting to pytest and widening option coverage are the
obvious next increments.

---

## 5. Packaging

- `conda_build/meta.yaml` was five releases stale at 3.23. Now 3.36, with the
  recipe test exercising both entry points' import graphs instead of one `-h`.
- `pyyaml`, `matplotlib-base` and `seaborn` are now declared. `pyyaml` is
  imported at module scope by `vsnp3_file_setup.py`, which both entry points
  import, so it was a hard requirement satisfied only incidentally by dask.
- A stale 1.4 GB Singularity image was committed. Bioconda builds from the GitHub
  source tarball, so every build and every user was downloading it. The tracked
  tree is now 4.8 MB. Publish container images as release assets.
- `conda_build/build.sh` now uses `set -euo pipefail` — a partly failed `mv`
  previously still reported a successful build — installs the shell helpers
  alongside the Python scripts, and sets the executable bit.
- `bin/get_stats.sh` and `bin/organize_directories.sh` were semicolon-chained, so
  an unguarded `cd ..` could be followed by `rm -r ./stats/` in the wrong
  directory, deleting a sibling run's data. Fixed before being installed.
- `internal/bump_version.sh` replaces the hand-edited `sed` recipe in
  `conda_build/instructions`.

### Dependency ranges, and the `<=` trap

`numpy` and `pandas` are no longer capped below the versions most of the
ecosystem has moved to. The old pins did not mean what they read as: **conda
pads the shorter version with zeros before comparing**, so `1.26.4 > 1.26` and
`numpy >=1.16,<=1.26` admitted numpy 1.26.0 and nothing else. The same applied
to `pandas <=2.2` (only 2.2.0) and `python >=3.8,<=3.12` (only 3.12.0).

That was not a harmless over-tightening. It held installed environments on numpy
1.26.0 while scipy 1.17 requires `>=1.26.4`, so `vsnp3_kernel_plots.py` was
calling `sns.kdeplot` — and through it `scipy.stats.gaussian_kde` — in a
configuration scipy declares unsupported. It also held users on Python 3.12.0,
from October 2023.

The new ranges are `numpy >=1.26.4,<3`, `pandas >=2.2.2,<3` and
`python >=3.12,<3.14`, with floors raised on `biopython`, `dask`, `sourmash`,
`matplotlib-base` and `seaborn`, and `scipy` newly declared — it was a hard
requirement of the kernel plots arriving only as a dependency of the `seaborn`
metapackage. `matplotlib-base` and `seaborn` had no floor at all, which left the
solver free to pair a numpy-1-ABI build with numpy 2.

Three things back the widening rather than assertion:

- No code in the repository uses a numpy API that numpy 2 changed. The one line
  that did — `int()` on a 1-D array in the N50 calculation — is fixed below.
- The seeded bootstrap draws identically on numpy 1.26.0, 2.2.6, 2.4.3 and
  2.4.6, so published support values do not move. `tests/test_numpy_contract.py`
  freezes that stream so a future numpy cannot change it quietly.
- CI now runs a version matrix at both ends of every declared range, and a
  `cross-version-agreement` job compares a full step2 run on numpy 1 against one
  on numpy 2, cell by cell, and fails on any changed call.

`python >=3.12` is a statement of fact rather than a preference:
`vsnp3_alignment_vcf.py` and `vsnp3_best_reference_sourmash.py` use PEP 701
f-strings, which are a `SyntaxError` before 3.12. The declared `>=3.8` would
have let conda install vSNP3 into an environment where step 1 dies on import.
`container/Singularity.def` pinned `python=3.11` and was in exactly that state.

`tests/test_release_consistency.py` now fails on any `<=` in the recipe, compares
`meta.yaml` and `conda_build/vsnp3.yaml` by version spec rather than by package
name, and requires a floor on every package that constrains the numpy ABI.

- **`bin/vsnp3_assembly.py`**: the N50 line called `int()` on the 1-D array
  returned by `np.where`. That is a `DeprecationWarning` up to numpy 2.2 and a
  hard `TypeError` from numpy 2.4 on, and it already raised on any assembly with
  a zero-length contig. The L50 line one below it had always indexed correctly.
  Step 1's assembly path is not covered by the end-to-end test, which is why it
  went unnoticed; `tests/test_numpy_contract.py` now covers it.

### Outstanding, for the maintainer

1. **`LICENSE` is missing.** `conda_build/meta.yaml` declares
   `license_file: LICENSE` and `license: GPL3`, but no such file has ever existed
   in this repository, and the badge in `README.md` links to one on GitHub. This
   is deliberately not auto-generated: for a US government work the choice
   between GPL3 and public domain under 17 USC §105 is a legal decision.
   `tests/test_release_consistency.py` will keep failing until it is added.
2. **Regenerate the `sha256`** in `conda_build/meta.yaml` against the 3.36 tag
   tarball before opening the bioconda PR. The environment has been confirmed to
   solve on `linux-64` and `osx-arm64` with the widened ranges, at the floor
   (numpy 1.26.4), at the current release, and alongside a pinned
   `numpy=2.1.3` — but re-check after any further dependency change.
3. **Notify collaborators running 3.35**, including ICAR-NIVEDI: every
   minus-strand amino acid they have produced is wrong.
4. **Fix the reference files.** `Brucella_ceti1` and `Brucella_suis4` have no gbk
   at all. `mtbc0_v1.1` has 4012 CDS with no `/product` and no `/translation`,
   which also means it gets no `self_test()` coverage. These are now reported at
   load time rather than failing quietly, but the files still need fixing.

### Deferred to 3.37

Sample-name truncation semantics, converting the remaining `os.system` calls in
`vsnp3_alignment_vcf.py`, `set -euo pipefail` across the `internal/` wrappers
(six of which currently contain Python docstrings pasted into bash and would
break immediately), per-sample scoping of step 1's temporary files, and an
opt-in `--strict` mode that turns the remaining loud fallbacks into errors.
