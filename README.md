🧬 Strain-Typer
Universal strain typing tool for Mash, FastANI, Skani, and any genome–genome comparison output

Strain-Typer is a flexible and universal tool that evaluates binning or clustering accuracy from any pairwise genome comparison table.
It supports:

Mash (distance → similarity auto-conversion)

FastANI

Skani

Any tabular comparison file (TSV)

Any reference mapping file containing isolate → strain cluster relationships

The tool auto-detects which type of score you provided (ANI, similarity, or distance) and makes accurate decisions accordingly.

🔧 Features

✔ Automatic score-type detection (ANI, similarity, distance)
✔ Supports any column order
✔ Detects isolate/substrain/strain columns in reference file
✔ Evaluates top-hit accuracy
✔ Produces:

Per-cluster accuracy plot

Cluster-to-cluster similarity heatmap

Summary report

Detailed output table
✔ Clean output folder per run
✔ Works with Mash / FastANI / Skani / sourmash / custom tools
