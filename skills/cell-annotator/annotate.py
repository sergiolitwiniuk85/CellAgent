"""
Cell Annotator — Marker-Based Cell Type Annotation
===================================================
Loads PanglaoDB marker database and annotates clusters by
overlap scoring between DEGs and known marker genes.
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional


def load_marker_db(
    db_path: str = "data/cell_markers/panglao_markers.tsv",
    species: str = "human",
    organ: Optional[str] = None,
    canonical_only: bool = False,
) -> Dict[str, List[str]]:
    """
    Load PanglaoDB markers into {cell_type: [genes]} dict.

    Parameters
    ----------
    db_path : str
        Path to PanglaoDB TSV file.
    species : str
        "human" → only human markers (Hs), "mouse" → only mouse (Mm),
        "all" → both species.
    organ : str or None
        Filter by organ (e.g., "Breast", "Immune system").
        None → all organs.
    canonical_only : bool
        If True, only include canonical markers (column 8 == 1).
    """
    path = Path(db_path)
    if not path.exists():
        raise FileNotFoundError(f"Marker database not found: {db_path}")

    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df.columns = df.columns.str.strip()

    # Filter by species
    if species == "human":
        df = df[df["species"].str.contains("Hs", na=False)]
    elif species == "mouse":
        df = df[df["species"].str.contains("Mm", na=False)]

    # Filter by organ
    if organ and organ.lower() != "all":
        df = df[df["organ"].str.lower() == organ.lower()]

    # Filter canonical only
    if canonical_only:
        df = df[df["canonical marker"] == "1"]

    # Build dict
    marker_dict = {}
    for _, row in df.iterrows():
        cell_type = str(row["cell type"]).strip()
        gene = str(row["official gene symbol"]).strip()
        if cell_type and gene and gene != "nan":
            marker_dict.setdefault(cell_type, []).append(gene)

    # Deduplicate genes per cell type
    return {ct: list(set(genes)) for ct, genes in marker_dict.items()}


def annotate_clusters(
    cluster_degs: Dict[str, List[str]],
    marker_dict: Dict[str, List[str]],
    min_overlap: int = 1,
    top_n: int = 50,
) -> Dict[str, Dict]:
    """
    Annotate clusters by overlapping DEGs with marker database.

    Parameters
    ----------
    cluster_degs : dict
        {cluster_id: [gene1, gene2, ...]}
    marker_dict : dict
        {cell_type: [gene1, gene2, ...]}
    min_overlap : int
        Minimum overlapping genes required to consider a match.
    top_n : int
        Number of top DEGs to use per cluster.

    Returns
    -------
    dict
        {cluster_id: {"cell_type": str, "score": float,
                       "matching_markers": [genes], "n_markers": int}}
    """
    results = {}
    for cluster, degs in cluster_degs.items():
        deg_set = set(degs[:top_n])
        best_type = "Unassigned"
        best_score = 0.0
        best_markers = []
        best_n = 0

        for cell_type, markers in marker_dict.items():
            marker_set = set(markers)
            overlap = deg_set & marker_set
            if len(overlap) >= min_overlap:
                score = len(overlap) / max(len(marker_set), 1)
                if score > best_score:
                    best_score = score
                    best_type = cell_type
                    best_markers = sorted(list(overlap))
                    best_n = len(overlap)

        results[cluster] = {
            "cell_type": best_type,
            "score": round(best_score, 3),
            "matching_markers": best_markers,
            "n_markers": best_n,
        }
    return results


def format_annotation_table(annotations: Dict[str, Dict]) -> str:
    """Format annotations as a markdown table."""
    lines = [
        "| Cluster | Predicted Cell Type | Score | Matching Markers |",
        "|---------|-------------------|-------|------------------|",
    ]
    for cluster in sorted(annotations.keys(), key=lambda x: int(x) if x.isdigit() else x):
        ann = annotations[cluster]
        markers_str = ", ".join(ann["matching_markers"][:5])
        if len(ann["matching_markers"]) > 5:
            markers_str += f" +{len(ann['matching_markers'])-5} more"
        if ann["score"] >= 0.5:
            score_str = f"**{ann['score']:.2f}** ★"
        elif ann["score"] >= 0.2:
            score_str = f"{ann['score']:.2f}"
        else:
            score_str = f"_{ann['score']:.2f}_"
        lines.append(
            f"| {cluster} | {ann['cell_type']} | {score_str} | {markers_str} |"
        )
    return "\n".join(lines)


def confidence_label(score: float) -> str:
    """Return confidence label for a score."""
    if score >= 0.5:
        return "High"
    elif score >= 0.2:
        return "Medium"
    elif score > 0:
        return "Low"
    return "No match"
