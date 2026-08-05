# Regression harness

Two properties of step 2 that are easy to break and impossible to notice by eye:

1. **Determinism** — the same input must produce byte-identical output. Before 3.36 it
   did not: the bootstrap RNG was unseeded, the annotated ALT at a multi-allele position
   was whichever VCF was read last, and several globs returned directory order. Two runs
   of the same six samples gave one clade 100% support and then 98%.
2. **Data identity across a refactor** — a performance change must not move a single
   call. A byte diff is the wrong tool for this, because legitimate ordering changes
   (sample order following filename order rather than directory order) rewrite every
   row and column without changing any value.

## Usage

```bash
export VSNP3_REGRESSION_DIR=/tmp/vsnp3_reg          # optional; a temp dir otherwise
DATA=~/vsnp3_test_dataset/AF2122_test_files/step2/original

tests/regression/run_step2.sh before "$DATA" Mycobacterium_AF2122 -html_tree
# ... make a change ...
tests/regression/run_step2.sh after  "$DATA" Mycobacterium_AF2122 -html_tree

tests/regression/compare.py "$VSNP3_REGRESSION_DIR"/snapshots/{before,after}
```

`compare.py` exits 0 when the data agrees and 1 when it does not, and separates the two
kinds of difference:

```
order  Mbovis-01__Mbovis-01.fasta: same records, different order
DATA   ..._position_tree.html.bootstrap: bootstrap support differs: ['100.0'] -> ['98.0']
```

`order` lines are presentation only. `DATA` lines are a real change and must be
explained before the change ships.

For determinism, run the same label twice and `diff -rq` the snapshots directly —
byte-identical is the bar, including the `.bootstrap` sidecars.

## What the tools do

- **`run_step2.sh`** copies the VCFs to a scratch directory before running. Not a
  convenience: step 2 deletes the files it reads once they are archived, and one of the
  things this harness checks is that it only does that to files it owns. Never point it
  at a pristine dataset.
- **`snapshot.py`** normalizes an output tree for comparison — strips run timestamps from
  filenames and bodies, converts every `.xlsx` to TSV, elides runtimes, and replaces
  plotly's per-render `uuid4` graph-div ids. Without the uuid rule every HTML tree looks
  different on every run and the real bootstrap variation is buried. Bootstrap
  percentages are also extracted to a `.bootstrap` sidecar so a diff can see them.
- **`compare.py`** compares SNP tables as `{(sheet, row, position) -> value}`, FASTA as
  `{sample -> sequence}`, newick as a leaf set, and bootstrap values as an ordered list.

Requires the reference set (see the README's Quick Start for
`vsnp3_test_dataset`) and `raxml` on `PATH`.
