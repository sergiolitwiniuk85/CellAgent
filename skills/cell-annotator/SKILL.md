---
name: cell-annotator
description: "Trigger: anotar, annotate, cell type, interpretación biológica, marker genes, lookup. Annotate cell clusters with cell-type labels using marker gene databases."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Cell Annotator — Marker-Based Cell Type Annotation

## Purpose

Assign biological cell-type labels to clusters using a marker gene lookup database. No LLMs, no GPUs, no paid APIs — pure overlap scoring between cluster DEGs and known marker genes from PanglaoDB.

## How It Works

```
Cluster DEGs (top 50 genes)
        ↓
    Match against PanglaoDB marker database (8,287 gene→cell-type entries)
        ↓
    For each cell type: compute overlap score = |DEGs ∩ Markers| / |Markers|
        ↓
    Assign best match with confidence = overlap / max_possible
        ↓
    Output: cluster → cell_type, score, matching_markers
```

## Database

PanglaoDB (March 2020 release):
- **8,287** gene→cell-type associations
- **259** unique cell types across **28** organs
- Columns: species, gene symbol, cell type, organ, sensitivity, specificity
- Stored at `data/cell_markers/panglao_markers.tsv`

## Annotation Logic

```python
def annotate_clusters(marker_dict: dict, cluster_degs: dict, min_overlap: int = 1):
    """
    marker_dict: {cell_type: [gene1, gene2, ...]}
    cluster_degs: {cluster_id: [gene1, gene2, ...]}
    
    Returns: {cluster_id: {"cell_type": str, "score": float, "markers": [genes]}}
    """
    results = {}
    for cluster, degs in cluster_degs.items():
        deg_set = set(deg[:50])  # top 50 DEGs
        best_type = "Unassigned"
        best_score = 0.0
        best_markers = []
        
        for cell_type, markers in marker_dict.items():
            marker_set = set(markers)
            overlap = deg_set & marker_set
            if len(overlap) >= min_overlap:
                # Jaccard-like score: overlap / len(markers)
                score = len(overlap) / max(len(marker_set), 1)
                if score > best_score:
                    best_score = score
                    best_type = cell_type
                    best_markers = list(overlap)
        
        results[cluster] = {
            "cell_type": best_type,
            "score": round(best_score, 3),
            "matching_markers": best_markers,
        }
    return results
```

## Usage

```python
from skills.cell_annotator.annotate import load_marker_db, annotate_clusters

# Load database
marker_db = load_marker_db("data/cell_markers/panglao_markers.tsv",
                           species="human", organ="all")

# Get DEGs per cluster (from scanpy)
degs = sc.get.rank_genes_groups_df(adata, group=None)
cluster_genes = {
    str(g): group["names"].tolist()
    for g, group in degs.groupby("group")
}

# Annotate
annotations = annotate_clusters(marker_db, cluster_genes)
```

## Output Format

```json
{
    "cluster_annotations": {
        "0": {"cell_type": "T cells", "score": 0.85, "matching_markers": ["CD3D", "CD7"]},
        "1": {"cell_type": "Macrophages", "score": 0.72, "matching_markers": ["CD68", "CD163"]},
        "2": {"cell_type": "Unassigned", "score": 0.0, "matching_markers": []}
    }
}
```

## Validation

- Score > 0.5 = High confidence
- Score 0.2-0.5 = Medium confidence (suggest review)
- Score < 0.2 = Low confidence (leave as "Unassigned")
- Clusters with 0 matching markers → report as "Unknown cluster type"

## Limitations

- Resolution limited to PanglaoDB's 259 cell types
- May not capture disease-specific subtypes (e.g., CAF subtypes like iCAF/myCAF)
- One-vs-many scoring — a gene can match multiple cell types
- Best for major lineages (immune, epithelial, stromal)

## Hard Rules

- Never fabricate annotations. If no markers match, say "Unassigned".
- Always include the confidence score alongside the label.
- Do not override the scientist's judgment — annotations are suggestions.
- The marker database is read-only. Never modify `data/cell_markers/`.
- **NEVER delete, overwrite, or modify the original input file**.
