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

📥 Installation
git clone https://github.com/yourusername/strain-typer.git
cd strain-typer
pip install -r requirements.txt

▶️ Usage
python3 strain_typer.py -i comparisons.tsv -r reference.tsv

Optional arguments:
--auto        # Automatically use columns 0,1,2 without prompts
--method XXX  # Tag the output folder (mash, fastani, etc)

📊 Example Outputs

After running, you will get:

results/method_timestamp/
│
├── strain_typer_output.tsv
├── run_summary.txt
├── accuracy_per_cluster.png
└── cluster_similarity_heatmap.png

📁 Input Format
Comparison file (TSV)

Can be Mash, FastANI, Skani, or custom format.

Must contain at least:

query   reference   score


Other columns are ignored.

Reference file (TSV)

Must contain columns representing:

isolate OR genome ID

substrain / sample ID

strain / strain cluster

The script will automatically detect them.

🧪 Example
python3 strain_typer.py \
  -i example_data/example_comparisons.tsv \
  -r example_data/example_reference.tsv \
  --method fastani
