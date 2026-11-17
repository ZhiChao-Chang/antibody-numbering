# Validation record — 2026-11-17

## Primary reference and denominators

The primary target is frozen abYsis software output. The representative cohort has 200 public inputs (100 H, 50 kappa, 50 lambda), of which 197 have reference numbering. The enriched difficult cohort has 100 public inputs, of which 95 have reference numbering. Reference-unclassified cases are tracked separately and excluded from the numbered-reference denominator.

| Scheme | Representative whole-input exact | Difficult whole-input exact |
| --- | ---: | ---: |
| Kabat | 197/197 | 91/95 |
| Chothia | 197/197 | 91/95 |
| IMGT | 197/197 | 91/95 |
| Martin | 197/197 | 91/95 |
| AHo | 197/197 | 89/95 |

Whole-input exactness includes assigned labels, blank positions, and domain coverage. Per-residue agreement uses reference-numbered residues as its denominator. Neither metric is experimental biological ground truth.

## Development relationship

These fixed sets have been opened. Public framework fragments and coordinate paths were curated from their SAbDab2/OAS references, and rules/thresholds were iterated against them. The table measures development compatibility, not untouched population accuracy. Inference receives sequence text and the bundled framework library; it never receives a per-query website answer.

## Remaining differences

Four difficult inputs have shorter returned domains in the website response. On the exact reference-returned domain, this implementation matches all 20 case/scheme groups; original whole-input metrics above remain unchanged. Two additional lambda first-loop cases differ only in AHo gap placement. The reference contains two layouts for the same loop length and Kabat path, and a justified general discriminator is not established. These cases are not classified as reference mistakes.

## Regression and expanded audit

- All 184 original regression records retain their mapping/rejection behavior. These are records, not 184 independent unique sequences; private input sequences are not distributed.
- Corrected negative controls: 1,000 records, zero classified/numbered outputs. Categories include ordinary proteins, non-antibody Ig-like proteins, TCR, and composition-matched shuffles.
- Expanded cohorts: diagnostic 1,000; locked 1,000; natural 10,000; synthetic boundary 500. All completed without runtime errors. Completion is separate from successful or complete numbering; per-scheme status counts are in the JSON summary.
- Combined audit: 13,984 records / 13,645 unique sequences after overlap and duplicates. Exact-score/search optimizations preserve every complete output relative to the pre-optimization candidate.
- All 300 primary website HTML receipts were re-parsed offline, with response hashes, label cardinality, unique labels, and domain offsets checked.

ANARCI/ANARCII cached outputs provide auxiliary implementation comparisons, not replacement truth. For the 10,000 natural inputs, post-iteration per-reference-residue agreement ranges from approximately 90.6% to 99.4% across the five schemes and two peers. Persistent convention differences remain: AHo whole-input exactness against both peers is 0 despite high residue agreement. These results are kept separate from the abYsis table.

## Reproduce the distributed checks

```bash
python3 -B -S -m unittest discover -s tests -v
python3 scripts/build_browser_manifest.py --check
```

The repository includes 33 unit, CLI, presentation/export, and public failure-mode regressions. The complete campaign inputs, raw reference responses, and private regression sequences are not bundled; the commands above reproduce the distributed checks, not the whole campaign. A future independent accuracy study needs a newly reserved cohort.

[validation_summary.json](validation_summary.json) contains aggregate results, scope, auxiliary comparisons, and [release.json](release.json) identifies exact published assets. The initial import inventory remains historical and is not rewritten.
