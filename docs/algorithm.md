# Algorithm and implementation guide

This guide explains the public implementation in execution order. The engine is a deterministic, reference-assisted sequence algorithm: it combines fixed search constraints, sequence alignment, and explicit numbering rules. It has no training phase or learned model weights.

## Where to start reading

| Component | Responsibility |
| --- | --- |
| `scripts/number_antibodies_independent.py::read_fasta()` | Parse and normalize FASTA while preserving original residue indices. |
| `engine.py::align()` | Compute global sequence alignment and residue correspondences. |
| `engine.py::AnchorNumberer.locate()` | Select a supported chain, framework template, and segmentation hypothesis. |
| `engine.py::AnchorNumberer.framework_alignment()` | Align individual framework segments and handle supported truncations. |
| `engine.py::loop_labels()` | Calculate scheme-specific loop labels and insertion positions. |
| `engine.py::AnchorNumberer.number()` | Assemble framework and loop mappings for all five schemes. |
| `scripts/number_antibodies_independent.py::FrequencyTables` | Optionally annotate computed positions with public amino-acid count frequencies. |

The engine is in [`src/antibody_numbering/engine.py`](../src/antibody_numbering/engine.py). Batch parsing, export, and provenance are in the [CLI script](../scripts/number_antibodies_independent.py).

## 1. Normalize input and retain coordinates

Each FASTA record represents one antibody variable domain. Parsing joins sequence lines and removes whitespace. Sequences are uppercased; alignment gaps (`-`) are removed and `?` becomes the unknown residue `X`. The input ledger retains the original sequence and a one-based mapping from normalized residues to original sequence offsets. Duplicate headers are preserved and distinguished using generated record IDs.

Records with stop characters, invalid residue symbols, or no sequence are skipped with an explicit status. Malformed FASTA structure raises an error. The engine's direct Python entry point performs residue validation but does not replace the CLI's complete parsing and provenance workflow.

## 2. Load public framework references

The default JSON asset contains five references: 6VJA VH/VL, 6XDG REGN10987 VL, 1N8Z VH, and 2VIR VL. Each stores framework sequences, zero-based half-open segment slices, chain/subtype metadata, and five sets of public framework labels. Loop sequences are masked with `X` and loop labels are blank.

These templates are development references, not independent validation examples. Input names and sequence hashes are not used to retrieve numbering answers. New loop labels are calculated by the rules below.

## 3. Locate a candidate framework

`locate()` evaluates candidate chain/template/segment combinations:

1. Candidate FR3 endpoints match at least two residues in a local `YYC` motif. This is a hard filter.
2. The first anchor is searched within the first 71 sequence offsets. Its expected `C` residue is a soft preference rather than a hard requirement.
3. Candidate segment lengths follow fixed supported ranges: the first inter-framework segment spans 10–12 residues for H or 10–17 for L; the second spans 15–19 for H or exactly 7 for L. Framework segment lengths impose additional constraints.
4. FR1, FR2, and FR3 are aligned separately. Candidates need at least 0.2 alignment-score units per residue in each segment and 0.45 across the framework segments together.
5. Conserved `C`/`W` residues contribute +6/+3. Recognized FR4 patterns add a +12/-12 preference for compatible/incompatible kappa or lambda templates. Missing pattern evidence contributes no subtype bonus.
6. The ranking score is the framework alignment total plus bonuses, divided by the number of framework residues. Ties prefer more framework residues, then follow deterministic chain/template ordering.

The lengths, motifs, scores, and thresholds are hand-set heuristics. They were not optimized through a training procedure and are not calibrated correctness probabilities. Failure to find a candidate returns `unclassified`; it does not establish that the sequence is biologically non-antibody.

## 4. Align and label framework residues

`align()` implements Needleman–Wunsch dynamic programming. The scoring model uses +3 for a match, -1 for a mismatch, 0 for pairs involving `X`, and -5 for insertions. Reference deletion penalties normally equal -5, with a supported FR1 terminal-truncation exception. Ties prefer diagonal, deletion, then insertion.

The traceback returns zero-based `(reference_index, query_index)` pairs, with `None` for a gap. Matched pairs transfer reference labels to query residues. Framework insertions without an implemented label rule raise `UnsupportedNumbering` for the affected scheme.

Heavy FR2 may use a different public framework sequence when its Kabat label signature is identical. Supported short FR4 tails are aligned from their start to avoid incorrectly shifting surviving residues to terminal labels.

## 5. Compute loop labels

The third inter-framework segment ends at the first supported FR4 motif. Without that motif, the remaining sequence is treated as an incomplete third segment. This is an algorithmic fallback, not evidence that biological FR4 is absent.

`loop_labels()` uses chain type, loop index, length, output scheme, terminal completeness, and light-chain subtype. Internal calculation spans can include flanking residues and must not be equated directly with every scheme's biological CDR definitions.

Alphabetic insertions occupy scheme-specific sites; for example, a supported first heavy segment can insert after Kabat H35 or Chothia/Martin H31. IMGT distributes shorter segments from both ends of an interval, leaving central gaps. Supported long third segments use decimal insertion strings around 111/112; labels are strings rather than floating-point values.

Explicit length checks and rule cases are in the source. An unsupported conversion returns a reason and blank positions for that scheme while preserving the other schemes' results.

### Legacy AHo compatibility

For a particular 19-residue H2 calculation segment, the verified legacy abYsis converter returns only the first six segment labels and truncates the remaining AHo mapping. The implementation preserves this behavior, leaves the suffix blank, and adds `compatibility_notes`. It is a compatibility rule and must not be interpreted as a complete biological numbering assignment.

## 6. Annotate frequencies and export results

Optional frequency tables are keyed by scheme and H/L-prefixed position. The percentage is `100 * amino_acid_count / total_including_gap`, using the public `All` organism category. Missing tables, unknown residues, and unmapped positions remain blank with explicit status. These lookups follow numbering and never influence framework selection or labels. No data are fetched during inference.

The CLI exports residue-level CSV, per-record/per-scheme summaries, the full inference result, an input ledger, missing frequency positions, and a run manifest. Code and reference hashes identify the exact assets used. Original input sequences remain in the ledger, so generated output must be handled with the same confidentiality as the input.

## Interpretation and validation boundaries

- Overall and per-scheme status must both be inspected. `numbered` indicates completion under the current rules and can coexist with blank positions or compatibility notes.
- The manifest's `OFFLINE_NUMBERING_FINISHED_REQUIRES_AUDIT` phase explicitly separates execution completion from independent validation of new inputs.
- The internal alignment score and optional reference frequency are neither confidence scores nor probabilities of correct numbering.
- Current support targets single variable domains within implemented framework and loop-length rules. Extreme loops, unusual insertions/deletions, divergent sequences, and multidomain constructs require additional validation and can be rejected or misassigned.
- The included tests check specific behavior and regressions. They do not establish population-level accuracy or universal equivalence with another numbering tool.
