# 🧬 Strain-Typer

**Strain-Typer** is a modular Python pipeline for strain-level classification
and benchmarking from pairwise genome comparison results.

It processes genome–genome comparison outputs from tools such as **FastANI,
Mash, and Skani**, integrates genome metadata and strain/substrain assignments,
identifies the best genomic match for each query genome, evaluates
classification accuracy, and generates automated visualizations.

> Strain-Typer is a post-processing and benchmarking tool. It does not replace
> the underlying genome-comparison algorithms.

## Features

- Process pairwise genome comparison results
- Supports FastANI, Mash, Skani, Sourmash-style, and custom tabular outputs
- Handles ANI, similarity, and distance score formats
- Automatic score-type detection with explicit override options
- Automatic detection of common query/reference/score columns
- Explicit column selection when needed
- Genome/isolate-to-strain metadata mapping
- Top-hit strain classification
- Overall and per-strain-cluster accuracy
- Mismatch reporting
- Per-cluster accuracy plot
- Cluster-to-cluster similarity heatmap
- Timestamped output directories
- Run summary and detailed output table

## Installation

```bash
git clone https://github.com/kmittal46/Strain_typer.git
cd Strain_typer
pip install -r requirements.txt
```

Python 3.9+ is recommended.

## Dependencies

- Python
- NumPy
- pandas
- Matplotlib

FastANI, Mash, Skani, and Sourmash do not need to be installed to process
their existing output files.

## Quick Start

A small example is included in `examples/`:

```bash
python strain_typer.py     -i examples/example_comparisons.tsv     -r examples/example_reference.tsv     --method fastani     --score-type ani
```

Results are written to:

```text
results/fastani_YYYYMMDD_HHMMSS/
├── strain_typer_output.tsv
├── run_summary.txt
├── accuracy_per_cluster.png
└── cluster_similarity_heatmap.png
```

## Input Files

### Pairwise comparison file

The comparison file should contain at least:

```text
query    reference    score
```

Additional columns are allowed.

Example:

```text
query    reference    ANI
genome_A genome_B     99.95
genome_A genome_C     99.91
```

Raw headerless Mash `mash dist` output is also supported.

### Reference metadata

The reference file connects genome/isolate identifiers to strain clusters:

```text
isolate    substrain    strain
genome_A   genome_A     cluster_1
genome_B   genome_B     cluster_1
genome_D   genome_D     cluster_2
```

Common names such as `isolate`, `genome`, `accession`, `substrain`,
`sample`, `strain`, and `cluster` can be detected automatically.

## Command-Line Usage

### Basic command

```bash
python strain_typer.py -i comparisons.tsv -r reference.tsv
```

- `-i` / `--input`: pairwise genome comparison file.
- `-r` / `--ref`: reference metadata file.

Without additional options, detected columns can be confirmed interactively.

### Automatic column selection

```bash
python strain_typer.py     -i comparisons.tsv     -r reference.tsv     --auto
```

`--auto` uses:

```text
column 0 = query
column 1 = reference
column 2 = score
```

This is useful for known, consistently formatted inputs.

### Specify the comparison method

```bash
python strain_typer.py     -i comparisons.tsv     -r reference.tsv     --method fastani
```

`--method` labels the output directory. It does not run the external comparison
program.

### Specify the score type

Use:

```text
--score-type ani
--score-type distance
--score-type similarity
--score-type auto
```

For ANI values such as `99.95`, `ani` converts the percentage to a 0–1
similarity:

```text
99.95 -> 0.9995
```

For distance scores:

```text
similarity = 1 - distance
```

For scores already expressed as similarity from 0 to 1, use `similarity`.

`auto` lets the program determine the score type.

### FastANI

```bash
python strain_typer.py     -i fastani_results.tsv     -r reference.tsv     --method fastani     --score-type ani
```

### Mash

Raw Mash `mash dist` output reports distance, where lower values indicate
greater similarity:

```bash
python strain_typer.py     -i mash_results.tsv     -r reference.tsv     --method mash     --score-type distance
```

Strain-Typer converts Mash distance using:

```text
similarity = 1 - distance
```

### Skani

For an ANI/similarity column:

```bash
python strain_typer.py     -i skani_results.tsv     -r reference.tsv     --method skani     --score-type ani
```

Choose the score type that matches the actual score column.

### Explicit column selection

For non-standard column names:

```bash
python strain_typer.py     -i comparisons.tsv     -r reference.tsv     --query-col query     --reference-col reference     --score-col ANI     --score-type ani
```

### Exclude self-matches

If query and reference sets contain the same genomes:

```bash
python strain_typer.py     -i comparisons.tsv     -r reference.tsv     --exclude-self
```

This removes comparisons where normalized query and reference IDs are identical.

For benchmarking, excluding self-matches is recommended when query and reference
sets overlap.

## Outputs

Each run creates:

```text
results/
└── method_YYYYMMDD_HHMMSS/
    ├── strain_typer_output.tsv
    ├── run_summary.txt
    ├── accuracy_per_cluster.png
    └── cluster_similarity_heatmap.png
```

### `strain_typer_output.tsv`

Contains the selected top reference hit for each query and the classification
information used by the pipeline.

### `run_summary.txt`

Records input information, detected columns, score information, mapping
statistics, classification statistics, and output locations.

### `accuracy_per_cluster.png`

Shows classification accuracy for each strain cluster. Bars are annotated with
correct/total predictions and percentage accuracy.

### `cluster_similarity_heatmap.png`

Shows mean normalized similarity between strain clusters. The visualization
focuses on the high-similarity range.

## Workflow

```text
Pairwise genome comparisons
            ↓
Input/column detection
            ↓
Score handling
            ↓
Normalize to similarity
            ↓
Genome → strain mapping
            ↓
Select top hit per query
            ↓
Predicted vs. query strain
            ↓
Accuracy metrics
            ↓
Tables + visualizations
```

See [`docs/workflow.md`](docs/workflow.md) for more detail.

## Example Dataset

The `examples/` directory contains a small synthetic dataset:

```text
examples/
├── example_comparisons.tsv
└── example_reference.tsv
```

Run it with:

```bash
python strain_typer.py     -i examples/example_comparisons.tsv     -r examples/example_reference.tsv     --method fastani     --score-type ani
```

The example contains multiple genomes assigned to the same strain clusters so
the workflow demonstrates actual strain-level classification.

## Benchmarking and Interpretation

Strain-Typer can be used to compare classification behavior across different
genome-comparison methods by processing their output files with the same
reference metadata.

Important considerations:

- Self-matches can artificially inflate accuracy.
- Query and reference genomes should ideally be independent.
- Each strain cluster should contain enough genomes for meaningful evaluation.
- Reference metadata determines the ground-truth strain assignment.
- Different comparison tools report different score semantics.
- Always verify the score type before interpreting results.

## Repository Structure

```text
Strain_typer/
├── strain_typer.py
├── README.md
├── requirements.txt
├── LICENSE
├── .gitignore
├── examples/
│   ├── example_comparisons.tsv
│   └── example_reference.tsv
├── results/
│   └── .gitkeep
└── docs/
    └── workflow.md
```

## Large-Scale Use

The repository intentionally does not contain large genome datasets or millions
of pairwise comparisons. Keep large input and result files outside GitHub and
use the included small example to demonstrate the workflow.

## License

This project is licensed under the MIT License. See `LICENSE` for details.
