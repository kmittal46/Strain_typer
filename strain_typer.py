#!/usr/bin/env python3
"""
Strain-Typer
============

A modular post-processing and benchmarking tool for pairwise genome
comparison results.

The script accepts comparison outputs from tools such as FastANI, Mash,
Skani, Sourmash, or custom tabular files; maps genome identifiers to
strain clusters; selects the best reference hit for each query; and
generates accuracy and similarity visualizations.

Example:
    python strain_typer.py \
        -i examples/example_comparisons.tsv \
        -r examples/example_reference.tsv \
        --method fastani \
        --score-type ani
"""

import argparse
import csv
import re
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# ID UTILITIES
# ============================================================

def normalize_id(value):
    """Normalize genome/isolate identifiers for metadata matching."""
    if pd.isna(value):
        return ""

    value = str(value).strip()
    value = value.split("/")[-1].split("\\")[-1]

    # Remove common FASTA extensions, including compressed variants.
    value = re.sub(
        r"\.(?:fa|fna|fasta)(?:\.gz)?$",
        "",
        value,
        flags=re.IGNORECASE,
    )
    return value


# ============================================================
# INPUT FILE HANDLING
# ============================================================

def detect_delimiter(path):
    """
    Detect the delimiter used by a text table.

    Returns one of:
        "\\t"  tab
        ","   comma
        None  whitespace
    """
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if "\t" in line:
                return "\t"
            if "," in line:
                return ","
            return None

    raise ValueError(f"Input file is empty: {path}")


def first_data_line(path):
    """Return the first non-empty, non-comment line."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                return line
    return ""


def split_line(line, delimiter):
    """Split a line according to the detected delimiter."""
    if delimiter == "\t":
        return line.split("\t")
    if delimiter == ",":
        return next(csv.reader([line]))
    return re.split(r"\s+", line.strip())


def looks_numeric(value):
    """Return True if a value can be interpreted as a number."""
    try:
        float(str(value).strip())
        return True
    except (TypeError, ValueError):
        return False


def looks_like_header(fields):
    """
    Detect common comparison-file headers.

    This prevents a headerless FastANI/Mash output from losing its first
    comparison row when pandas reads it.
    """
    known_terms = {
        "query",
        "reference",
        "query_name",
        "reference_name",
        "ani",
        "identity",
        "similarity",
        "score",
        "distance",
        "dist",
        "mashdist",
        "fragments_mapped",
        "fragments_total",
        "p-value",
        "pvalue",
        "shared-hashes",
        "shared_hashes",
    }

    lowered = {str(field).strip().lower() for field in fields}

    if lowered & known_terms:
        return True

    # A conventional header should generally contain at least one
    # non-numeric field in a position where a score is expected.
    if len(fields) >= 3 and not looks_numeric(fields[2]):
        return True

    return False


def identify_headerless_format(fields):
    """
    Identify common headerless genome-comparison formats.

    FastANI:
        query reference ANI fragments_mapped fragments_total

    Mash dist:
        query reference distance p_value shared_hashes
    """
    if len(fields) >= 5:
        # Mash's fifth field is typically something like 10000/10000.
        if re.fullmatch(r"\d+/\d+", str(fields[4]).strip()):
            return "mash"

        # FastANI's fourth/fifth fields are numeric fragment counts.
        if looks_numeric(fields[2]) and looks_numeric(fields[3]) and looks_numeric(fields[4]):
            return "fastani"

    if len(fields) >= 3 and looks_numeric(fields[2]):
        return "generic"

    return None


def load_comparison_file(path):
    """
    Load a comparison table while preserving headerless raw FastANI/Mash
    outputs.

    Returns:
        dataframe, detected_format, delimiter
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Comparison file not found: {path}")

    delimiter = detect_delimiter(path)
    line = first_data_line(path)

    if not line:
        raise ValueError(f"Comparison file is empty: {path}")

    fields = split_line(line, delimiter)
    has_header = looks_like_header(fields)

    if has_header:
        if delimiter == "\t":
            df = pd.read_csv(path, sep="\t", dtype=str)
        elif delimiter == ",":
            df = pd.read_csv(path, sep=",", dtype=str)
        else:
            df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)

        detected_format = "headered"
    else:
        detected_format = identify_headerless_format(fields)

        if delimiter == "\t":
            sep = "\t"
        elif delimiter == ",":
            sep = ","
        else:
            sep = r"\s+"

        df = pd.read_csv(
            path,
            sep=sep,
            engine="python",
            header=None,
            dtype=str,
        )

        if detected_format == "mash":
            names = [
                "Query",
                "Reference",
                "Distance",
                "PValue",
                "SharedHashes",
            ]
            if df.shape[1] > len(names):
                names.extend(
                    [f"Column_{i}" for i in range(len(names), df.shape[1])]
                )
            df.columns = names[: df.shape[1]]

        elif detected_format == "fastani":
            names = [
                "Query",
                "Reference",
                "ANI",
                "Fragments_Mapped",
                "Fragments_Total",
            ]
            if df.shape[1] > len(names):
                names.extend(
                    [f"Column_{i}" for i in range(len(names), df.shape[1])]
                )
            df.columns = names[: df.shape[1]]

        else:
            df.columns = [
                f"Column_{i}" for i in range(df.shape[1])
            ]

    return df, detected_format, delimiter


# ============================================================
# SCORE HANDLING
# ============================================================

def detect_score_type(series, column_name="", method="run"):
    """
    Detect one of:
        ani
        distance
        similarity
        unknown

    Column names and method names are preferred over numeric heuristics.
    Explicit --score-type should be used whenever the score semantics are
    ambiguous.
    """
    column_lower = str(column_name).lower()
    method_lower = str(method).lower()

    # Strong column-name signals.
    if "ani" in column_lower or "identity" in column_lower:
        return "ani"

    if any(term in column_lower for term in ("distance", "dist", "mashdist")):
        return "distance"

    if "similarity" in column_lower:
        return "similarity"

    # Strong method-specific signal for raw Mash distance output.
    if method_lower == "mash":
        return "distance"

    numeric = pd.to_numeric(series, errors="coerce").dropna()

    if len(numeric) == 0:
        return "unknown"

    smin = float(numeric.min())
    smax = float(numeric.max())
    median = float(numeric.median())

    # ANI values are normally reported as percentages.
    if 70 < smax <= 100:
        return "ani"

    # For values already in [0, 1], a high-valued distribution is most
    # consistent with similarity while a low-valued distribution is most
    # consistent with distance.
    if 0 <= smin and smax <= 1:
        if median >= 0.5:
            return "similarity"
        return "distance"

    return "unknown"


def normalize_scores(series, score_type):
    """Convert ANI, distance, or similarity to a 0-1 similarity scale."""
    numeric = pd.to_numeric(series, errors="coerce")

    if score_type == "ani":
        similarity = numeric / 100.0
    elif score_type == "distance":
        similarity = 1.0 - numeric
    elif score_type == "similarity":
        similarity = numeric
    else:
        similarity = numeric

    return similarity.clip(lower=0.0, upper=1.0)


# ============================================================
# COLUMN DETECTION
# ============================================================

def find_named_column(columns, exact_names, partial_names):
    """Find a column using exact matches first, then partial matches."""
    columns = list(columns)
    lowered = [str(column).strip().lower() for column in columns]

    for name in exact_names:
        if name in lowered:
            return columns[lowered.index(name)]

    for name in partial_names:
        for index, column_name in enumerate(lowered):
            if name in column_name:
                return columns[index]

    return None


def detect_comparison_columns(df):
    """Detect query, reference, and score columns."""
    columns = list(df.columns)

    query_col = find_named_column(
        columns,
        exact_names=[
            "query",
            "query_name",
            "query_name_clean",
            "qname",
            "queryid",
            "query_id",
            "q_id",
            "q",
        ],
        partial_names=["query", "qname", "q_id"],
    )

    reference_col = find_named_column(
        columns,
        exact_names=[
            "reference",
            "ref",
            "ref_name",
            "reference_name",
            "refname",
            "rname",
            "referenceid",
            "reference_id",
            "r",
        ],
        partial_names=["reference", "ref", "rname"],
    )

    score_col = find_named_column(
        columns,
        exact_names=[
            "ani",
            "identity",
            "similarity",
            "score",
            "distance",
            "dist",
            "mashdist",
        ],
        partial_names=["ani", "identity", "similarity", "distance", "dist", "score"],
    )

    return query_col, reference_col, score_col


def choose_column(columns, detected, label, default_index):
    """Interactively confirm a detected column or request an index."""
    if detected is not None:
        use = input(
            f"Detected column '{detected}' that looks like {label}. "
            "Use it? [Y/n]: "
        ).strip().lower()

        if use in ("", "y", "yes"):
            return detected

    try:
        index = int(
            input(
                f"Index for {label.upper()} column "
                f"(default {default_index}): "
            )
            or default_index
        )
        return columns[index]
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Invalid {label} column selection.") from exc


# ============================================================
# REFERENCE METADATA
# ============================================================

def detect_ref_columns(ref_df):
    """
    Detect isolate, substrain, and strain cluster columns.

    Exact names are preferred, followed by common alternatives and finally
    positional fallbacks for backwards compatibility.
    """
    columns = list(ref_df.columns)

    iso_col = find_named_column(
        columns,
        exact_names=["isolate"],
        partial_names=["isolate", "genome", "accession", "iso"],
    )

    sub_col = find_named_column(
        columns,
        exact_names=["substrain", "sub_strain"],
        partial_names=["substrain", "sub_strain", "sample"],
    )

    strain_col = find_named_column(
        columns,
        exact_names=["strain", "strain_cluster", "cluster"],
        partial_names=[
            "strain_cluster",
            "strain_representative",
            "strain_rep",
            "strain",
            "cluster",
        ],
    )

    # Preserve the original positional fallback behavior.
    if iso_col is None and len(columns) >= 1:
        iso_col = columns[0]
    if sub_col is None and len(columns) >= 2:
        sub_col = columns[1]
    if strain_col is None and len(columns) >= 3:
        strain_col = columns[2]

    return iso_col, sub_col, strain_col


def load_reference_file(path):
    """Load and validate the reference metadata table."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Reference file not found: {path}")

    delimiter = detect_delimiter(path)

    if delimiter == "\t":
        ref_df = pd.read_csv(path, sep="\t", dtype=str)
    elif delimiter == ",":
        ref_df = pd.read_csv(path, sep=",", dtype=str)
    else:
        ref_df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)

    if ref_df.empty:
        raise ValueError(f"Reference file is empty: {path}")

    return ref_df


def build_strain_map(ref_df, iso_col, sub_col, strain_col):
    """Build identifier -> strain-cluster mapping."""
    ref_map = {}

    for _, row in ref_df.iterrows():
        strain = (
            str(row[strain_col]).strip()
            if pd.notna(row[strain_col])
            else ""
        )

        iso = (
            normalize_id(row[iso_col])
            if pd.notna(row[iso_col])
            else ""
        )
        sub = (
            normalize_id(row[sub_col])
            if pd.notna(row[sub_col])
            else ""
        )

        if iso:
            ref_map[iso] = strain

        if sub:
            ref_map[sub] = strain

    return ref_map


# ============================================================
# CLASSIFICATION
# ============================================================

def get_top_hits(valid, exclude_self=False):
    """
    Select the highest-similarity reference for each query.

    Similarity is already normalized so higher values are always better.
    """
    working = valid.copy()

    if exclude_self:
        working = working[
            working["Genome2"].str.lower()
            != working["Genome1"].str.lower()
        ].copy()

    top = (
        working
        .sort_values(
            ["Genome2", "Similarity"],
            ascending=[True, False],
            kind="mergesort",
        )
        .drop_duplicates("Genome2")
        .copy()
    )

    top["Correct"] = top["Strain1"] == top["Strain2"]

    return top


def calculate_accuracy(top):
    """Return total classifications, correct classifications, and accuracy."""
    total = len(top)
    correct = int(top["Correct"].sum()) if total else 0
    accuracy = (correct / total * 100) if total else 0.0
    return total, correct, accuracy


# ============================================================
# VISUALIZATION
# ============================================================

def plot_cluster_accuracy(top, outdir, log):
    """Generate the per-strain-cluster accuracy plot."""
    log("Generating improved accuracy plot...")

    if top.empty:
        log("No top-hit classifications available; skipping accuracy plot.")
        return

    accuracy = (
        top.groupby("Strain2")["Correct"]
        .mean()
        .sort_values(ascending=False)
        * 100
    )

    counts = (
        top.groupby("Strain2")
        .size()
        .reindex(accuracy.index)
        .fillna(0)
        .astype(int)
    )

    correct_counts = (
        top.groupby("Strain2")["Correct"]
        .sum()
        .reindex(accuracy.index)
        .fillna(0)
        .astype(int)
    )

    # Preserve the original threshold-based colors.
    colors = []
    for value in accuracy:
        if value < 90:
            colors.append("#D62728")
        elif value < 95:
            colors.append("#FF8C00")
        else:
            colors.append("#2CA02C")

    n_clusters = len(accuracy)

    # Scale the figure for larger numbers of clusters.
    width = min(max(8, n_clusters * 0.25), 80)

    fig, ax = plt.subplots(figsize=(width, 6))
    bars = ax.bar(
        range(n_clusters),
        accuracy.values,
        color=colors,
        edgecolor="black",
    )

    # Keep annotations for small/medium plots. For very large cluster
    # counts, labels become unreadable, so retain the bars without overlap.
    annotate = n_clusters <= 60

    if annotate:
        for index, bar in enumerate(bars):
            height = bar.get_height()
            correct = int(correct_counts.iloc[index])
            total = int(counts.iloc[index])
            percentage = accuracy.iloc[index]

            label = f"{correct}/{total}\n{percentage:.0f}%"
            y_position = min(height + 1.2, 99.5)

            ax.text(
                bar.get_x() + bar.get_width() / 2,
                y_position,
                label,
                ha="center",
                va="bottom",
                fontsize=8 if n_clusters <= 30 else 6,
            )

    ax.set_ylim(0, 102)
    ax.set_xticks(range(n_clusters))
    ax.set_xticklabels(
        accuracy.index,
        rotation=90,
        fontsize=8 if n_clusters <= 60 else 6,
    )
    ax.set_ylabel("Accuracy (%)", fontsize=12)
    ax.set_xlabel("Strain cluster", fontsize=12)
    ax.set_title("Accuracy Per Strain Cluster", fontsize=14)

    fig.tight_layout()
    output_path = outdir / "accuracy_per_cluster.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    log(f"Saved {output_path.name}")


def plot_similarity_heatmap(valid, outdir, log):
    """Generate the cluster-to-cluster mean similarity heatmap."""
    log("Generating cluster similarity heatmap...")

    if valid.empty:
        log("No valid comparisons available; skipping heatmap.")
        return

    similarity_matrix = valid.pivot_table(
        index="Strain2",
        columns="Strain1",
        values="Similarity",
        aggfunc="mean",
    )

    similarity_matrix = similarity_matrix.reindex(
        index=sorted(similarity_matrix.index),
        columns=sorted(similarity_matrix.columns),
    )

    sim_values = similarity_matrix.to_numpy(dtype=float)

    nrows, ncols = sim_values.shape

    # Scale for larger matrices while preventing enormous output images.
    fig_height = min(max(6, nrows * 0.12), 120)
    fig_width = min(max(6, ncols * 0.12), 120)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Keep missing cells visually distinct instead of treating missing
    # comparisons as zero similarity.
    cmap = plt.colormaps["magma"].copy()
    cmap.set_bad("#E6E6E6")

    masked_values = np.ma.masked_invalid(sim_values)

    image = ax.imshow(
        masked_values,
        cmap=cmap,
        vmin=0.90,
        vmax=1.00,
        aspect="auto",
    )

    colorbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.04,
        pad=0.02,
    )
    colorbar.set_label("Mean Similarity", fontsize=12)

    ax.set_xticks(range(ncols))
    ax.set_xticklabels(
        similarity_matrix.columns,
        rotation=90,
        fontsize=6,
    )

    ax.set_yticks(range(nrows))
    ax.set_yticklabels(
        similarity_matrix.index,
        fontsize=6,
    )

    ax.set_title(
        "Cluster-to-Cluster Mean Similarity Heatmap",
        fontsize=14,
    )

    fig.tight_layout()

    output_path = outdir / "cluster_similarity_heatmap.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    log(f"Saved {output_path.name}")


# ============================================================
# OUTPUT / LOGGING
# ============================================================

def create_output_directory(method):
    """Create a timestamped output directory."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path("results") / f"{method}_{timestamp}"
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def make_logger(summary_file):
    """Create a logger that prints to screen and writes to the summary."""
    def log(message=""):
        print(message)
        with open(summary_file, "a", encoding="utf-8") as handle:
            handle.write(str(message) + "\n")

    return log


def write_top_hits(top, outdir):
    """Write the top-hit classification table."""
    output_path = outdir / "strain_typer_output.tsv"
    top.to_csv(output_path, sep="\t", index=False)
    return output_path


# ============================================================
# COMMAND-LINE INTERFACE
# ============================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Classify genomes by their highest-similarity reference "
            "and evaluate strain-level accuracy."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Pairwise genome comparison file.",
    )

    parser.add_argument(
        "-r",
        "--ref",
        required=True,
        help=(
            "Reference metadata file mapping isolate/substrain "
            "identifiers to strain clusters."
        ),
    )

    parser.add_argument(
        "--auto",
        action="store_true",
        help=(
            "Use columns 0, 1, and 2 as query, reference, and score "
            "without interactive prompts."
        ),
    )

    parser.add_argument(
        "--query-col",
        help="Name of the query column.",
    )

    parser.add_argument(
        "--reference-col",
        help="Name of the reference column.",
    )

    parser.add_argument(
        "--score-col",
        help="Name of the score column.",
    )

    parser.add_argument(
        "--score-type",
        choices=["auto", "ani", "distance", "similarity"],
        default="auto",
        help=(
            "Interpretation of the comparison score. "
            "Default: auto."
        ),
    )

    parser.add_argument(
        "--method",
        default="run",
        help=(
            "Comparison method name used for labeling the output "
            "directory (e.g. fastani, mash, skani)."
        ),
    )

    parser.add_argument(
        "--exclude-self",
        action="store_true",
        help="Exclude comparisons where query and reference IDs are identical.",
    )

    return parser.parse_args()


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    args = parse_args()

    # --------------------------------------------------------
    # Output setup
    # --------------------------------------------------------
    outdir = create_output_directory(args.method)
    summary_file = outdir / "run_summary.txt"
    log = make_logger(summary_file)

    log(f"Input: {args.input}")
    log(f"Ref:   {args.ref}")
    log("")

    # --------------------------------------------------------
    # Load comparison data
    # --------------------------------------------------------
    df, detected_format, delimiter = load_comparison_file(args.input)

    log(f"Detected input format: {detected_format}")
    log(
        "Detected delimiter: "
        + (
            "TAB"
            if delimiter == "\t"
            else "COMMA"
            if delimiter == ","
            else "WHITESPACE"
        )
    )
    log(f"Loaded {len(df):,} rows.")

    columns = df.columns.tolist()

    log("Columns detected:")
    for index, column in enumerate(columns):
        log(f"[{index}] {column}")

    # --------------------------------------------------------
    # Select comparison columns
    # --------------------------------------------------------
    detected_q, detected_r, detected_score = detect_comparison_columns(df)

    if args.query_col:
        if args.query_col not in columns:
            raise ValueError(
                f"Query column '{args.query_col}' was not found. "
                f"Available columns: {columns}"
            )
        q_col = args.query_col
    elif args.auto:
        if len(columns) < 3:
            raise ValueError(
                "--auto requires at least three columns "
                "(query, reference, score)."
            )
        q_col = columns[0]
    else:
        q_col = choose_column(columns, detected_q, "query", 0)

    if args.reference_col:
        if args.reference_col not in columns:
            raise ValueError(
                f"Reference column '{args.reference_col}' was not found. "
                f"Available columns: {columns}"
            )
        r_col = args.reference_col
    elif args.auto:
        r_col = columns[1]
    else:
        r_col = choose_column(columns, detected_r, "reference", 1)

    if args.score_col:
        if args.score_col not in columns:
            raise ValueError(
                f"Score column '{args.score_col}' was not found. "
                f"Available columns: {columns}"
            )
        s_col = args.score_col
    elif args.auto:
        s_col = columns[2]
    else:
        s_col = choose_column(columns, detected_score, "score", 2)

    log(
        f"\nUsing -> QUERY: {q_col} ; "
        f"REFERENCE: {r_col} ; SCORE: {s_col}"
    )

    # --------------------------------------------------------
    # Normalize IDs and scores
    # --------------------------------------------------------
    df["Genome2"] = df[q_col].apply(normalize_id)
    df["Genome1"] = df[r_col].apply(normalize_id)

    df["__raw_score__"] = df[s_col].astype(str)
    df["Score_numeric"] = pd.to_numeric(
        df[s_col],
        errors="coerce",
    )

    if args.score_type == "auto":
        score_type = detect_score_type(
            df["Score_numeric"],
            column_name=s_col,
            method=args.method,
        )
    else:
        score_type = args.score_type

    log("\nSample raw Score values (first 10):")
    log(df["__raw_score__"].head(10).tolist())

    non_numeric = int(df["Score_numeric"].isna().sum())

    log(
        f"Non-numeric after coercion: "
        f"{non_numeric} / {len(df):,}"
    )

    log(f"\nDetected score type: {score_type}")

    if score_type == "unknown":
        raise ValueError(
            "Could not determine the score type automatically. "
            "Please use --score-type ani, distance, or similarity."
        )

    df["Similarity"] = normalize_scores(
        df["Score_numeric"],
        score_type,
    )

    # --------------------------------------------------------
    # Score statistics
    # --------------------------------------------------------
    log("\nScore stats (post-coercion to numeric):")
    log(str(df["Score_numeric"].describe()))

    if df["Similarity"].notna().any():
        log(
            "Similarity range: "
            f"{df['Similarity'].min():.6f} .. "
            f"{df['Similarity'].max():.6f}"
        )
    else:
        log("Similarity range: no valid numeric values")

    # --------------------------------------------------------
    # Load reference metadata
    # --------------------------------------------------------
    ref_df = load_reference_file(args.ref)

    iso_col, sub_col, strain_col = detect_ref_columns(ref_df)

    if not iso_col or not sub_col or not strain_col:
        raise ValueError(
            "Could not identify isolate, substrain, and strain columns "
            "in the reference file."
        )

    log(
        f"\nReference columns detected: "
        f"isolate='{iso_col}', "
        f"substrain='{sub_col}', "
        f"strain='{strain_col}'"
    )

    ref_map = build_strain_map(
        ref_df,
        iso_col,
        sub_col,
        strain_col,
    )

    # --------------------------------------------------------
    # Apply reference mapping
    # --------------------------------------------------------
    df["Strain1"] = df["Genome1"].map(ref_map)
    df["Strain2"] = df["Genome2"].map(ref_map)

    mapped1 = int(df["Strain1"].notna().sum())
    mapped2 = int(df["Strain2"].notna().sum())

    valid = df.dropna(
        subset=["Strain1", "Strain2", "Similarity"]
    ).copy()

    log(
        f"\nMapped genome1 -> strain: {mapped1:,} ; "
        f"genome2 -> strain: {mapped2:,}"
    )

    log(
        "Valid mapped comparisons "
        f"(both mapped & similarity present): {len(valid):,}"
    )

    if args.exclude_self:
        before = len(valid)
        valid = valid[
            valid["Genome2"].str.lower()
            != valid["Genome1"].str.lower()
        ].copy()

        log(
            f"Self-comparisons excluded: "
            f"{before - len(valid):,}"
        )

    log("\nSample mappings (first 10 valid rows):")
    for _, row in valid.head(10).iterrows():
        log(
            f"{row['Genome2']} -> {row['Strain2']}   "
            f"(ref {row['Genome1']} -> {row['Strain1']})  "
            f"sim={row['Similarity']:.4f}"
        )

    # --------------------------------------------------------
    # Top-hit classification
    # --------------------------------------------------------
    top = get_top_hits(
        valid,
        exclude_self=False,
    )

    total, correct, accuracy = calculate_accuracy(top)

    log("\n--- Top-hit summary ---")
    log(f"Top-hit comparisons: {total:,}")
    log(f"Correct: {correct:,}")
    log(f"Accuracy: {accuracy:.2f}%")

    output_path = write_top_hits(top, outdir)
    log(f"Saved top-hits to: {output_path}")

    mismatches = top[top["Correct"] == False]

    if len(mismatches) > 0:
        log("\nSample mismatches (first 10):")
        for _, row in mismatches.head(10).iterrows():
            log(
                f"  Query:{row['Genome2']} -> "
                f"Ref:{row['Genome1']}  "
                f"sim={row['Similarity']:.4f}  "
                f"cluster_query={row['Strain2']} "
                f"cluster_ref={row['Strain1']}"
            )

    # --------------------------------------------------------
    # Visualizations
    # --------------------------------------------------------
    plot_cluster_accuracy(top, outdir, log)
    plot_similarity_heatmap(valid, outdir, log)

    # --------------------------------------------------------
    # Finish
    # --------------------------------------------------------
    log("\nDone.")
    log(f"Results folder: {outdir}")


if __name__ == "__main__":
    main()
