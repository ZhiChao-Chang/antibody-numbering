# Public reference data

This directory contains the framework asset required for offline numbering and optional snapshots of public position-frequency tables. Both are local inference inputs; the tool does not submit sequences to a website.

## Framework asset

`independent_numbering_public_20261109/framework_templates.json` uses the `public_numbered_frameworks_v2` schema and contains five templates:

| Template | Public source | Chain |
| --- | --- | --- |
| `H` | [6VJA](https://www.rcsb.org/structure/6VJA), rituximab | VH |
| `L` | [6VJA](https://www.rcsb.org/structure/6VJA), rituximab | Kappa VL |
| `L_lambda` | [6XDG](https://www.rcsb.org/structure/6XDG), REGN10987 | Lambda VL |
| `H_1N8Z` | [1N8Z](https://www.rcsb.org/structure/1N8Z), Herceptin | VH |
| `L_lambda_mouse` | [2VIR](https://www.rcsb.org/structure/2VIR) | Mouse lambda VL |

Each template provides a sequence, chain/subtype metadata, framework segment boundaries, and labels for Kabat, Chothia, IMGT, Martin, and AHo. Segment boundaries are zero-based half-open slices. Framework labels provide positional landmarks. Loop sequences are masked with `X` and loop labels are blank; the engine calculates the new input's loop labels.

These are development references, not independent holdout controls. They neither require a new antibody to have the same target nor restrict numbering to antibodies from those complexes. Source statements remain embedded in the JSON.

## Expanded library (2026-11-17)

The same v2 asset additionally contains `framework_profiles` (152 public framework fragments) and `framework_position_paths` (38 five-scheme coordinate paths). Fragments are FR1/FR2/FR3/FR4 pieces, not complete antibodies or CDR answer sequences. Each retains chain/locus, labels, source, and public source ID.

The library was curated from opened SAbDab2/OAS development examples with cached abYsis annotations. Selection checks annotation consistency, covers new label topologies, and uses 0.86 positional identity for sequence coverage; runtime scoring uses independent alignment thresholds. See [framework_library_provenance.json](framework_library_provenance.json), [algorithm.md](../../docs/algorithm.md), and [validation.md](../../docs/validation.md). Original base templates and optional frequency snapshots are retained.

## Frequency snapshots

`abysis_frequencies_all_20261105/` and `independent_numbering_frequencies_20261109/` contain JSON tables named by scheme and public position, for example `kabat_H1.json`. The source is [abYsis](http://www.abysis.org/abysis/). Each table retains:

- Numbering scheme, H/L-prefixed position, and organism category (`All`).
- Per-residue counts and the total including gap observations.
- Source URL, retrieval timestamp, and original response SHA-256.
- The stored website display percentages.

The CLI calculates percentages from counts rather than copying rounded display text. The snapshots cover a subset of positions and may not cover a new input. A missing reference is left blank and reported explicitly. Files in the first supplied frequency directory take precedence; no table is automatically refreshed or downloaded.

Reference frequency describes observed amino-acid usage at a position. It does not estimate the probability that an input's numbering assignment is correct.

## Provenance and rights

The [initial import inventory](../../docs/provenance.json) records the imported files' SHA-256 hashes. A new run separately records the framework and frequency references actually used in `manifest.json`.

Public structures, numbering labels, and frequency materials retain their original sources and applicable rights. The repository's MIT code license does not relicense third-party reference materials. No partner sequences or private experiment outputs are included in this directory.
