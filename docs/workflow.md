# Strain-Typer Workflow

Strain-Typer is a post-processing and benchmarking layer for pairwise genome
comparison results. Genome comparisons are performed by external tools such as
FastANI, Mash, or Skani; Strain-Typer integrates those results with strain metadata.

## Workflow

```text
Genome comparison tool
        |
        v
Pairwise comparison output
        |
        v
   Strain-Typer
        |
        +--> Input and column detection
        |
        +--> Score normalization
        |      ANI -> 0-1 similarity
        |      distance -> 1-distance
        |      similarity -> unchanged
        |
        +--> Genome -> strain mapping
        |
        +--> Top reference hit per query
        |
        +--> Predicted vs. query strain
        |
        +--> Accuracy calculation
        |
        +--> Visualizations + summary
```

## Input

Comparison files should contain a query genome, reference genome, and score.
Additional columns are allowed.

Reference metadata maps genome/isolate identifiers to strain clusters.

## Classification

For each query, valid comparisons are ranked by normalized similarity. The
highest-similarity reference is selected as the top hit. The prediction is
correct when the reference strain matches the query strain in the metadata.

## Self-matches

If query and reference sets overlap, a genome can match itself with a perfect
score. For benchmarking, use `--exclude-self` when appropriate.

## Outputs

Each run creates a timestamped directory under `results/` containing:

- `strain_typer_output.tsv`
- `run_summary.txt`
- `accuracy_per_cluster.png`
- `cluster_similarity_heatmap.png`

Generated results are ignored by Git by default.

