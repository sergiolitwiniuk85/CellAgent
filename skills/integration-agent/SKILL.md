---
name: integration-agent
description: "Trigger: integración, multimodal, muon, WNN, MOFA, CITE-seq, ATAC+RNA, multi-omics. Run multimodal integration with muon."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Integration Agent — Integración Multimodal

## Propósito

Integrar datos multimodales (RNA + ATAC, RNA + proteína, etc.) usando muon y herramientas del ecosistema scverse.

## Formato de Entrada

El agente espera un `MuData` con dos o más modadalidades:

```python
mdata = mu.read("path/to/data.h5mu")
# mdata.mod["rna"]   → AnnData con RNA
# mdata.mod["atac"]  → AnnData con ATAC (o "prot" para proteína)
```

## Pipeline de Integración

### 1. Preprocesamiento por Modalidad

**RNA**: usar normalize-agent (si no se hizo antes)

```python
sc.pp.normalize_total(mdata["rna"], target_sum=1e4)
sc.pp.log1p(mdata["rna"])
sc.pp.highly_variable_genes(mdata["rna"], n_top_genes=2000, flavor="seurat_v3")
sc.pp.pca(mdata["rna"], n_comps=30, svd_solver="arpack")
```

**Proteína (CITE-Seq)**:

```python
import muon as mu

# Normalización específica para proteínas: CLR (centered log-ratio)
mu.prot.pp.clr(mdata["prot"])
# También se puede usar el módulo directo:
# from muon import prot as mp
# mp.pp.clr(mdata["prot"])

sc.pp.pca(mdata["prot"], n_comps=15, svd_solver="arpack")
```

**ATAC**:

```python
from muon import atac as ac

# LSI (Latent Semantic Indexing) — equivalente a PCA para ATAC
ac.tl.lsi(mdata["atac"], n_comps=30)
```

### 2. Weighted Nearest Neighbor (WNN)

WNN integra modadalidades pesando la contribución de cada una por célula:

```python
import muon as mu

# WNN analysis — la joya de muon
mu.tl.wnn(mdata, modality_weights="uniform")  # o pasar weights manual

# Clustering sobre el grafo WNN
sc.tl.leiden(mdata, neighbors_key="wnn", key_added="leiden_wnn", resolution=0.8)
```

### 3. MOFA+ (Multi-Omics Factor Analysis)

Alternativa a WNN: factores latentes compartidos entre modadalidades:

```python
# MOFA+ requiere instalación aparte: pip install mofapy2
# O usar desde R: MOFA2

# Nota: MOFA+ es más lento pero modela explícitamente la estructura
# de cada modadalidad (gaussiana para RNA, Bernoulli para ATAC, etc.)
```

### 4. UMAP conjunto

```python
# Sobre representación WNN
sc.tl.umap(mdata, neighbors_key="wnn", random_state=42)
sc.pl.umap(mdata, color=["leiden_wnn"], legend_loc="on data")

# O sobre factores MOFA
# sc.tl.umap(mdata, use_rep="X_MOFA")
```

### 5. Visualización Conjunta

```python
# Mostrar contribución de cada modadalidad
mu.pl.wnn(mdata)

# UMAP coloreado por modadalidad y clusters
sc.pl.umap(mdata, color=["leiden_wnn", "rna:n_genes_by_counts", "prot:total_counts"],
           legend_loc="on data")
```

## Interpretación para el Usuario

- **WNN**: asigna pesos por célula — si una célula tiene mejor señal de RNA, el RNA pesa más
- **MOFA+**: encuentra "factores" latentes compartidos entre modadalidades
- **Leiden en WNN**: clusters que usan info de TODAS las modadalidades

## Output

```python
{
    "stage": "integration",
    "status": "completed",
    "summary": f"Integración WNN completada. {n_clusters} clusters conjuntos",
    "method": "WNN",
    "modalities": ["rna", "prot"],
    "n_clusters": int(mdata.obs["leiden_wnn"].nunique()),
    "data_path": "/path/to/integrated_data.h5mu",
    "plots": ["wnn_weights.png", "umap_wnn_clusters.png", "umap_by_modality.png"],
    "recommendations": [
        "WNN muestra buena separación de tipos celulares",
        "¿Querés extraer marker genes por cluster multimodal?",
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
