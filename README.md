<h1 align="center">Antibody Numbering</h1>

<p align="center">
  <strong>Antibody variable-domain numbering in pure Python.</strong><br>
  FASTA sequences to five numbering schemes, in your browser or from the command line.
</p>

<p align="center">
  <a href="#quick-start"><img src="https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&amp;logo=python&amp;logoColor=white" alt="Python 3.10 or later"></a>
  <a href="#how-it-works"><img src="https://img.shields.io/badge/Runtime-CPU%20%7C%20Offline-167D7F?style=flat-square" alt="CPU and offline runtime"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/Code%20license-MIT-243B53?style=flat-square" alt="Code license: MIT"></a>
</p>

<p align="center">
  <a href="https://zhichao-chang.github.io/antibody-numbering/">Online tool</a> ·
  <a href="#quick-start">Quick start</a> ·
  <a href="#usage">Usage</a> ·
  <a href="#outputs">Outputs</a>
</p>

---

Assign **Kabat, Chothia, IMGT, Martin, and AHo** positions to antibody VH and VL sequences. Process single or multi-record FASTA, inspect residue-level mappings and CDR regions, and download results as CSV or ZIP.

## Online tool

**[Open Antibody Numbering →](https://zhichao-chang.github.io/antibody-numbering/)**

Paste your FASTA, click **Start numbering**, then choose a sequence and numbering scheme. CSV contains the full batch; ZIP organizes results into one directory per sequence. CDR definitions can be selected independently of the numbering scheme.

Numbering runs in your browser. Sequences stay in page memory until you clear the input or close the tab. The first run downloads approximately 13 MB of runtime files.

## Quick start

Requires **Python 3.10+**. Run directly from a checkout using the Python standard library.

```bash
git clone https://github.com/ZhiChao-Chang/antibody-numbering.git
cd antibody-numbering

python3 scripts/number_antibodies_independent.py \
  --input examples/6vja_fv.fasta \
  --output results/example
```

The example contains public rituximab VH/VL sequences from [PDB 6VJA](https://www.rcsb.org/structure/6VJA). Open `numbering_with_frequencies.csv` in the output directory to view the mappings. Use a new output directory for each run.

## Usage

Use **one VH or VL variable domain per FASTA record**. Heavy and light chains can share a file as separate records.

```bash
python3 scripts/number_antibodies_independent.py \
  --input antibodies.fasta \
  --output results/my_run
```

| Option | Description |
| --- | --- |
| `--input PATH` | Input FASTA file. |
| `--output PATH` | New directory for results. |
| `--frameworks PATH` | Optional custom framework JSON; the bundled library is the default. |
| `--frequency-cache PATH` | Optional local amino-acid frequency tables; may be repeated. |
| `--help` | Show command-line help. |

<details>
<summary>Optional amino-acid reference frequencies</summary>

```bash
python3 scripts/number_antibodies_independent.py \
  --input antibodies.fasta \
  --output results/with_frequencies \
  --frequency-cache data/reference/abysis_frequencies_all_20261105 \
  --frequency-cache data/reference/independent_numbering_frequencies_20261109
```

`frequency_percent` describes observed amino-acid usage at a numbered position. Tables cover a subset of positions; unavailable values remain blank. The browser tool uses numbering without frequency annotation.

</details>

## Outputs

| File | Contents |
| --- | --- |
| `numbering_with_frequencies.csv` | All sequences and five schemes, one row per residue and scheme. |
| `records_summary.csv` | Status and mapped-residue count for each sequence and scheme. |
| `backend_outputs.json` | Complete mappings, selected framework, segments, and compatibility notes. |
| `input_ledger.json` | Original/normalized sequences and their index correspondence. |
| `manifest.json` | Run status, record counts, input and software/reference hashes. |
| `missing_frequency_positions.json` | Positions without optional frequency tables. |

Indices are **1-based**. Keep standard labels such as `H35A` and `H112.1` as strings.

`unclassified` indicates an unrecognized domain, `unsupported` an unsupported scheme case, and `unmapped_position` a residue without a label. Invalid input records are skipped with a reason. Partial results have blank positions and compatibility notes, including the supported legacy AHo partial-conversion cases.

## How it works

Conserved sequence landmarks guide framework alignment. A bundled reference library and scheme-specific rules produce the five mappings. Browser and command-line versions use the same Python engine. Intended inputs are single antibody variable domains; unusual scaffolds, extreme loops, and multidomain constructs may require separate handling.

## Development

```bash
python3 -B -S -m unittest discover -s tests -v
python3 scripts/build_browser_manifest.py --check
```

See [web/README.md](web/README.md) for local preview and deployment. Report issues with a minimal public or synthetic example and the command or browser used.

## Data and license

Bundled framework references include public PDB, SAbDab2, and OAS sources with abYsis numbering annotations. See [data/reference/README.md](data/reference/README.md) for details.

Original code is licensed under [MIT](LICENSE). Third-party reference data and the bundled browser runtime retain their source attribution and licenses.
