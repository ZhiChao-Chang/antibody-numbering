<h1 align="center">Antibody Numbering</h1>

<p align="center">
  <strong>Offline antibody variable-domain numbering in pure Python.</strong><br>
  From FASTA sequences to five standard numbering schemes, with residue-level mappings and traceable outputs.
</p>

<p align="center">
  <a href="#quick-start"><img src="https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&amp;logo=python&amp;logoColor=white" alt="Python 3.10 or later"></a>
  <a href="#how-it-works"><img src="https://img.shields.io/badge/Runtime-CPU%20%7C%20Offline-167D7F?style=flat-square" alt="CPU and offline runtime"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/Code%20license-MIT-243B53?style=flat-square" alt="Code license: MIT"></a>
</p>

<p align="center">
  <a href="#quick-start">Quick start</a> ·
  <a href="#usage">Usage</a> ·
  <a href="#outputs">Outputs</a> ·
  <a href="#how-it-works">Method</a> ·
  <a href="#validation-and-scope">Validation</a>
</p>

---

Antibody Numbering assigns **Kabat, Chothia, IMGT, Martin, and AHo** positions to antibody variable-domain sequences. The implementation combines conserved-anchor detection, framework alignment, and explicit loop-numbering rules. It runs locally using the Python standard library and bundled public framework references.

| Capability | What you get |
| :--- | :--- |
| **Batch FASTA input** | Process heavy and light variable domains in a single file. |
| **Five numbering schemes** | Generate all five mappings in one invocation, including supported insertion labels. |
| **Residue traceability** | Retain input identifiers, original indices, normalized indices, and processing status. |
| **Offline execution** | Run on a CPU with no third-party Python packages, model downloads, or sequence submissions. |
| **Optional reference frequencies** | Annotate positions using bundled snapshots of public amino-acid counts. |

## Quick start

**Requirements:** Git and Python 3.10 or later. Run directly from a source checkout; no package installation is needed.

```bash
git clone https://github.com/ZhiChao-Chang/antibody-numbering.git
cd antibody-numbering

python3 -S scripts/number_antibodies_independent.py \
  --input examples/6vja_fv.fasta \
  --output results/example_001
```

The example contains the public rituximab VH and VL sequences from [PDB 6VJA](https://www.rcsb.org/structure/6VJA). It produces **10 sequence–scheme summaries** and **1,135 residue rows**. Open `results/example_001/numbering_with_frequencies.csv` for the mappings and `records_summary.csv` for per-sequence status.

Choose a **new output directory** for each run. The `-S` option disables third-party site-package loading; the tool only uses the standard library.

## Usage

### Number your sequences

```bash
python3 -S scripts/number_antibodies_independent.py \
  --input /path/to/antibodies.fasta \
  --output results/my_run
```

Use one antibody variable domain per FASTA record. Put VH and VL in separate records; both may be present in the same file. Record names are preserved and do not determine chain classification.

For complete antibody chains, extract the variable domain before numbering. General target proteins and concatenated multidomain constructs are outside the intended input scope.

### Command-line options

| Option | Required | Description |
| :--- | :---: | :--- |
| `--input PATH` | Yes | FASTA file containing one or more variable-domain sequences. |
| `--output PATH` | Yes | Output directory to create; must not already exist. |
| `--frameworks PATH` | No | Override the bundled framework JSON with a compatible reference asset. |
| `--frequency-cache PATH` | No | Add a local frequency-table directory. Repeat to search multiple directories in order. |
| `--help` | No | Display the command-line reference. |

The default framework file is resolved relative to the project, so it also works when the script is called by absolute path from another working directory:

```text
data/reference/independent_numbering_public_20261109/framework_templates.json
```

### Add amino-acid reference frequencies

```bash
python3 -S scripts/number_antibodies_independent.py \
  --input examples/6vja_fv.fasta \
  --output results/example_with_frequencies \
  --frequency-cache data/reference/abysis_frequencies_all_20261105 \
  --frequency-cache data/reference/independent_numbering_frequencies_20261109
```

`frequency_percent` is the observed amino-acid count divided by the reference total, including gaps, multiplied by 100. The bundled tables use the **All** organism category. They are fixed snapshots with partial position coverage; missing references and unknown residues receive a blank frequency and an explicit status. No tables are fetched during inference.

**Reference frequency is not a probability that a numbering assignment is correct.** Omitting the cache options leaves frequencies blank and still produces numbering results.

## Outputs

Each run writes six files:

| File | Purpose |
| :--- | :--- |
| `numbering_with_frequencies.csv` | Residue-level mappings for every sequence and scheme, with optional reference frequencies. |
| `records_summary.csv` | Processing status and mapped-residue count for each sequence–scheme pair. |
| `backend_outputs.json` | Full mappings, framework selection, inferred segments, and compatibility notes. |
| `input_ledger.json` | Original sequences, normalized sequences, and source-index mappings. |
| `manifest.json` | Input, source-code, and reference hashes; run statistics and frequency provenance. |
| `missing_frequency_positions.json` | Positions without a matching frequency reference. |

Selected columns from the public example:

```csv
header,scheme,chain_type,raw_index,input_residue,position,status
public_6VJA_H,kabat,H,1,Q,H1,numbered
public_6VJA_H,kabat,H,2,V,H2,numbered
public_6VJA_H,kabat,H,3,Q,H3,numbered
```

`raw_index` and `normalized_index` are **1-based**. Preserve `position` as a string: labels can include prefixes, letters, and decimal insertions, such as `H35A` and `H112.1`. The full CSV also includes the original header, normalized residue, frequency status, and reference path.

<details>
<summary><strong>Input normalization and status interpretation</strong></summary>

- Sequences are normalized to uppercase. Alignment gaps (`-`) are removed, and `?` is represented as the unknown residue `X`; original indices remain available.
- Records containing `*`, invalid characters, or an empty sequence are skipped with a reason in the ledger and summary.
- `unclassified` means no acceptable framework interpretation was found.
- `unsupported` means a scheme encountered a case beyond the implemented rules; details are recorded in `backend_outputs.json`.
- `unmapped_position` identifies residues without an assigned standard position.
- `numbered` records completion under the implemented rules. It does not certify correctness or guarantee that every residue received a label.

Malformed FASTA structure, such as sequence text before the first header, is rejected. Check the per-scheme results and compatibility notes when interpreting partial outputs.

</details>

## How it works

```mermaid
flowchart LR
    A[FASTA input] --> B[Normalize and locate anchors]
    B --> C[Select and align frameworks]
    C --> D[Map framework positions]
    C --> E[Compute loop labels]
    D --> F[Five numbering schemes]
    E --> F
    F --> G[CSV and JSON outputs]
```

1. **Locate the variable-domain framework.** Candidate segments are evaluated using conserved-residue evidence, length constraints, and framework similarity. The selected interpretation provides the heavy/light assignment and a reference framework.
2. **Align framework residues.** A locally implemented Needleman–Wunsch algorithm maps the input framework segments to numbered reference positions, with explicit handling for supported terminal truncations.
3. **Compute loop labels.** Scheme-specific rules use chain type, light-chain subtype, loop length, and terminal completeness to assign loop positions and supported insertion labels.
4. **Export the mappings.** All five schemes are written with processing status and provenance. Optional frequency annotation is applied after numbering and does not influence the assignments.

The five public templates comprise two VH and three VL frameworks from **6VJA, 6XDG, 1N8Z, and 2VIR**. Their loop sequences are masked and their loop labels are removed. These are development references; they are not an independent validation set.

Inference uses no ANARCI/ANARCII imports, HMMER calls, neural weights, or per-query website numbering responses. No training step is required. The internal alignment score is a selection score, not a calibrated confidence estimate.

For a function-by-function walkthrough, scoring details, coordinate conventions, and compatibility behavior, see the [implementation guide](docs/algorithm.md) and the annotated [core engine](src/antibody_numbering/engine.py).

## Python integration

The core engine can also be called from a Python process. Run this example from the repository root:

```python
import sys
from pathlib import Path

sys.path.insert(0, str(Path("src").resolve()))
from antibody_numbering import AnchorNumberer

numberer = AnchorNumberer(
    Path("data/reference/independent_numbering_public_20261109/framework_templates.json")
)

# Public rituximab VH sequence from PDB 6VJA.
sequence = (
    "QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAIYPGNGDTSYN"
    "QKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWGAGTTVTVSA"
)
result = numberer.number(sequence)

print(result["chain"])                       # H
print(result["schemes"]["imgt"]["positions"])  # One label per input residue; unmapped entries are blank.
```

The engine returns dictionaries of per-scheme results. For arbitrary inputs, inspect the overall and per-scheme status before accessing mappings; unclassified results have a different shape. Use the CLI for batch FASTA parsing, input normalization records, and frequency annotation.

## Validation and scope

Run the bundled test suite without additional dependencies:

```bash
python3 -B -S -m unittest discover -s tests -v
```

The current suite contains **18 tests** covering alignment, insertion labels, kappa/lambda-specific rules, terminal truncation, unknown residues, negative controls, offline inference, frequency arithmetic, default framework resolution, custom reference overrides, and output-directory protection.

**Current scope:** single antibody variable domains within the implemented framework and loop-length rules. Extreme loop lengths, unusual framework insertions or deletions, divergent sequences, and multidomain constructs may be rejected or assigned incorrectly. The bundled tests provide regression coverage, not a population-level accuracy estimate.

A legacy abYsis AHo conversion behavior for a particular 19-residue H2 segment is intentionally reproduced: only part of the mapping is returned. Such cases are recorded in `compatibility_notes`. Review these notes before using partial mappings in downstream analyses.

## Repository layout

```text
antibody-numbering/
├── scripts/number_antibodies_independent.py   # Batch CLI
├── src/antibody_numbering/                   # Core engine
├── data/reference/                          # Frameworks and optional frequency tables
├── examples/6vja_fv.fasta                    # Public VH/VL example
├── tests/                                   # Unit and CLI tests
└── docs/provenance.json                     # SHA-256 inventory of the initial import
```

Run manifests describe the assets used in each invocation. The provenance inventory records the initial imported bytes and is not a live checksum manifest for subsequent source edits. Generated results and local input directories are excluded by `.gitignore`.

## Contributing

[Open an issue](https://github.com/ZhiChao-Chang/antibody-numbering/issues) with the command, Python version, commit identifier, observed behavior, and expected result. Use a minimal public or synthetic sequence when sharing a reproducer. For numbering-rule changes, include the scheme, a traceable reference assignment, and a regression test.

## Reference data and license

Framework sources: [6VJA](https://www.rcsb.org/structure/6VJA), [6XDG](https://www.rcsb.org/structure/6XDG), [1N8Z](https://www.rcsb.org/structure/1N8Z), and [2VIR](https://www.rcsb.org/structure/2VIR). Optional position-frequency snapshots originate from [abYsis](http://www.abysis.org/abysis/); each table retains its source URL, retrieval timestamp, counts, and response hash.

The [reference-data guide](data/reference/README.md) documents the framework schema, table fields, and source provenance.

The original source code is available under the [MIT License](LICENSE). Public structures and third-party reference data retain their original attribution and applicable rights; the code license does not relicense those materials.
