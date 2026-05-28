---
name: integration-agent
description: "Trigger: integración, multimodal, muon, WNN, MOFA, CITE-seq, ATAC+RNA, multi-omics. Run multimodal integration with muon."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Integration Agent — Multimodal Integration

## Purpose

Integrate multimodal data (RNA + ATAC, RNA + protein, etc.) using muon and the scverse ecosystem.

## Input Format

The agent expects a `MuData` with two or more modalities:

```python
mdata = mu.read("path/to/data.h5mu")
# mdata.mod["rna"]   → AnnData with RNA
# mdata.mod["atac"]  → AnnData with ATAC (or "prot" for protein)
```

## Integration Pipeline

### 1. Per-Modality Preprocessing

**RNA**: use normalize-agent (if not already done)

```python
sc.pp.normalize_total(mdata["rna"], target_sum=1e4)
sc.pp.log1p(mdata["rna"])
sc.pp.highly_variable_genes(mdata["rna"], n_top_genes=2000, flavor="seurat_v3")
sc.pp.pca(mdata["rna"], n_comps=30, svd_solver="arpack")
```

**Protein (CITE-Seq)**:

```python
import muon as mu

# Protein-specific normalization: CLR (centered log-ratio)
mu.prot.pp.clr(mdata["prot"])
# Alternative direct module:
# from muon import prot as mp
# mp.pp.clr(mdata["prot"])

sc.pp.pca(mdata["prot"], n_comps=15, svd_solver="arpack")
```

**ATAC**:

```python
from muon import atac as ac

# LSI (Latent Semantic Indexing) — PCA equivalent for ATAC
ac.tl.lsi(mdata["atac"], n_comps=30)
```

### 2. Weighted Nearest Neighbor (WNN)

WNN integrates modalities by weighting each modality's contribution per cell:

```python
import muon as mu

# WNN analysis — the crown jewel of muon
mu.tl.wnn(mdata, modality_weights="uniform")  # or pass manual weights

# Clustering on the WNN graph
sc.tl.leiden(mdata, neighbors_key="wnn", key_added="leiden_wnn", resolution=0.8)
```

### 3. MOFA+ (Multi-Omics Factor Analysis)

Alternative to WNN: shared latent factors across modalities:

```python
# MOFA+ requires separate install: pip install mofapy2
# Or from R: MOFA2

# Note: MOFA+ is slower but explicitly models each modality's
# data structure (Gaussian for RNA, Bernoulli for ATAC, etc.)
```

### 4. Joint UMAP

```python
# On WNN representation
sc.tl.umap(mdata, neighbors_key="wnn", random_state=42)
sc.pl.umap(mdata, color=["leiden_wnn"], legend_loc="on data")

# Or on MOFA factors
# sc.tl.umap(mdata, use_rep="X_MOFA")
```

### 5. Joint Visualization

```python
# Show contribution of each modality
mu.pl.wnn(mdata)

# UMAP colored by modality and clusters
sc.pl.umap(mdata, color=["leiden_wnn", "rna:n_genes_by_counts", "prot:total_counts"],
           legend_loc="on data")
```

## Interpretation for the User

- **WNN**: assigns per-cell weights — if a cell has better RNA signal, RNA weighs more
- **MOFA+**: finds shared latent "factors" across modalities
- **Leiden on WNN**: clusters using info from ALL modalities

## Output

```python
{
    "stage": "integration",
    "status": "completed",
    "summary": f"WNN integration complete. {n_clusters} joint clusters",
    "method": "WNN",
    "modalities": ["rna", "prot"],
    "n_clusters": int(mdata.obs["leiden_wnn"].nunique()),
    "data_path": "/path/to/integrated_data.h5mu",
    "plots": ["wnn_weights.png", "umap_wnn_clusters.png", "umap_by_modality.png"],
    "recommendations": [
        "WNN shows good cell type separation",
        "Would you like to extract marker genes per multimodal cluster?",
    ]
}
```

## Hard Rules

- Each modality needs its own preprocessing BEFORE integration
- Do not mix scales: RNA is log-normalized, protein is CLR, ATAC is LSI
- Save as .h5mu (not .h5ad) to preserve the multimodal structure
- If a modality has fewer than 100 cells, warn that WNN may fail
- **NEVER delete, overwrite, or modify the original input file**. Integrated data is always saved to a NEW output path.
- The only person who removes original data is the human from the CLI.
