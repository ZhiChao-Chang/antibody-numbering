# Reference data

## Framework library

`independent_numbering_public_20261109/framework_templates.json` is the default numbering reference. It contains framework sequences, chain metadata, and coordinate mappings for the five numbering schemes. It is loaded automatically by both the CLI and browser tool.

Sources include public PDB structures [6VJA](https://www.rcsb.org/structure/6VJA), [6XDG](https://www.rcsb.org/structure/6XDG), [1N8Z](https://www.rcsb.org/structure/1N8Z), and [2VIR](https://www.rcsb.org/structure/2VIR), plus SAbDab2/OAS framework fragments with abYsis numbering annotations. Individual source IDs are retained in the JSON.

## Optional frequency tables

`abysis_frequencies_all_20261105/` and `independent_numbering_frequencies_20261109/` contain [abYsis](http://www.abysis.org/abysis/) position-frequency snapshots, named by scheme and position, such as `kabat_H1.json`.

Pass one or both directories with `--frequency-cache` to annotate observed amino-acid usage. Percentages use residue counts divided by the total including gaps, for organism category `All`. Missing positions remain blank. Each table retains its source URL, retrieval date, counts, and response hash.

## Attribution

Public structures and third-party reference materials retain their original attribution and applicable rights. The repository's MIT license covers its original code.
