# Algorithm and implementation

This is a reference-assisted, deterministic variable-domain numberer implemented with the Python standard library. It targets compatibility with the saved abYsis annotations. Runtime accepts sequence text and bundled reference assets; external numbering tools and website responses are not queried during inference.

## 1. Input and coordinate conventions

The CLI reads one VH/VL variable domain per FASTA record. It removes whitespace and `-`, uppercases residues, and maps `?` to unknown `X`. Stops, invalid symbols, and empty records receive explicit skip statuses. Duplicate headers remain separate records. An index map preserves correspondence with the original sequence.

Engine offsets are zero-based, half-open slices. CSV raw and normalized indices are one-based. Labels such as `H35A`, `H112.1`, and `L106A` are strings, not numerical array indices or decimal values.

## 2. Reference library

The default v2 JSON asset contains:

- Five masked base frameworks from 6VJA, 6XDG, 1N8Z, and 2VIR.
- 152 public framework fragments: FR1 80, FR2 24, FR3 32, FR4 16.
- 38 coordinate paths keyed by chain, light subtype, framework region, and its ordered Kabat labels.

Fragments were curated from public SAbDab2/OAS inputs with cached abYsis labels. Selection groups chain/locus and five-scheme label topology, prioritizes new coordinate layouts, and greedily covers sequence classes using 0.86 direct positional identity. Identical fragments and conflicting annotations are checked before selection. This is reference-library curation and development calibration, not neural training.

The runtime fragment bank contains framework sequence and labels, not CDR sequences or complete antibody-answer records. Base-template loops are masked with `X` and their labels are blank. Sources and selection settings are retained in the asset and [library provenance](../data/reference/framework_library_provenance.json).

## 3. Candidate segmentation

`AnchorNumberer.locate()` enumerates chain/template/anchor/length hypotheses. The first anchor is searched within the first 71 offsets. A C contributes support but is not mandatory. Candidate FR3 endpoints meet at least two of three local features: Y/F, Y/F, C.

| Search span | Heavy | Light |
| --- | --- | --- |
| First internal loop | 7–16 | 8–20 |
| Second internal loop | 15–24 | 7–15 |
| FR2 | 13 or 14 | 15 or 16 |
| FR3 | 27–36 | 30–40 |

These are internal calculation spans, which can include flanking residues; they do not define every scheme's biological CDR boundaries.

FR1, FR2, and FR3 each require >=0.2 score units per residue; their combined unbonused mean must be >=0.45. Matches score +3, mismatches -1, unknown X pairs 0, and ordinary gaps -5. C and W anchor support adds +6 and +3. Recognized compatible light-chain terminal patterns add +12; contradictory kappa/lambda types are rejected.

The ranking score is `(sum(framework_scores) + anchor_bonus + subtype_bonus) / framework_residues`, minus `0.05 * max(0, default_FR2_length - candidate_FR2_length)` and `0.03 * max(0, default_FR3_length - candidate_FR3_length)`. Defaults are H 14/30 and L 15/32. Candidates need a final score >=1.3. Ties favor fuller frameworks and a deterministic ordering. These hand-set, development-calibrated heuristics are not probabilities.

## 4. Framework matching

`align()` implements Needleman–Wunsch dynamic programming with explicit traceback pairs and deterministic diagonal/deletion/insertion tie handling. `ordered_framework_alignment()` directly maps equal-length segments and uses deletion-only dynamic programming for shorter queries; longer queries retain general insertion evidence. Supported N-terminal truncation uses reduced leading deletion penalties (-0.1); short FR4 tails keep prefix coordinates.

`framework_choice()` considers matching chain/locus fragments. A valid donor must improve the base score by at least 4 for FR1/FR4 or 8 for FR2/FR3. FR2 donors must match the explicitly searched segment length. An unhandled query insertion cannot silently receive a label.

## 5. FR4 and domain completion

`select_fr4()` evaluates multiple terminal windows using full-prefix sequence evidence, rather than selecting the first motif inside CDR3. Candidates need mean score >=1.2; tails shorter than seven residues require an F/W start. A recognizable legacy motif is retained as a fallback for divergent tails. The selected boundary ends the third internal loop. If no boundary is found, the remaining tail follows incomplete-segment rules.

## 6. Framework coordinates

Most FR1/FR2/FR4 labels follow the selected alignment, then use a known framework coordinate path when available. A deleted Kabat slot need not imply the same label gap in IMGT/AHo. For example, one 22-residue light FR1 path omits Kabat L10 but IMGT L7.

`framework3_labels()` constructs chain/length-specific FR3 coordinates for heavy chains and light spans of at least 32 residues. It preserves terminal cysteine anchors and places extras at scheme-specific sites, e.g. heavy Kabat H82 versus Martin H72. Shorter light FR3 spans retain alignment-specific handling.

## 7. Loop labels and compatibility

`loop_labels()` uses chain, subtype, loop index/length, scheme, and terminal completeness. A 12-residue first heavy segment uses H35A/H35B in Kabat but H31A/H31B in Chothia/Martin. IMGT places short loops from both interval ends and supports explicit insertion strings for longer segments. Kappa/lambda AHo gap layouts have separate cases.

Unsupported conversions return `unsupported` with blank positions and a reason for that scheme; other schemes are preserved. A supported second heavy segment of length >=19 intentionally reproduces observed legacy abYsis AHo truncation: only six H2-segment labels survive, the following suffix is blank, and `compatibility_notes` records the behavior.

## 8. Exact computational optimizations

Conserved first-anchor candidates are evaluated early. An optimistic bound of +3 per remaining framework residue, plus maximum possible bonuses, prunes candidates unable to beat an existing result; original tie behavior is preserved.

Equal-length framework fragments use five-bit residue fields, an X mask, XOR, and integer `bit_count()` to reproduce exactly the +3/-1/0 substitution sum. Traceback pairs are allocated only for winning donors. This changes computation, not the scoring model. Randomized and exhaustive single-residue tests cover equivalence, and complete outputs match the pre-optimization candidate on 13,645 unique inputs.

## 9. Exports and frequency annotation

The CLI emits a residue CSV, per-record/scheme summaries, complete engine results, an input ledger, missing-frequency records, and source/input/asset hashes. `numbered` can include unmapped positions or compatibility notes; inspect the per-scheme data. The browser uses this same CLI and then applies display/export helpers.

Optional reference frequencies are looked up after numbering as `100 * residue_count / total_including_gap`. Missing values remain blank with explicit statuses. Frequencies cannot affect numbering and are not correctness probabilities.

## 10. Validation interpretation

The public compatibility cohorts were opened during development and informed the framework library and rules. Four difficult differences concern reference-returned domain coverage; all 20 case/scheme groups match when supplied exactly that returned domain. Two remaining AHo first-loop gap layouts have no established general selection rule. They remain unresolved. See [validation.md](validation.md) for denominators and engineering checks.
