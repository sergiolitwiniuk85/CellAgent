---
name: normalize-agent
description: "Trigger: normalizar, normalización, HVG, highly variable genes, scaling, log1p, normalize. Run normalization and feature selection for single-cell data."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Normalize Agent — Normalización y Selección de Features

## Propósito

Ejecutar la normalización estándar de scanpy: conteos por célula, log-transform, selección de genes altamente variables (HVG), y escalado.

## Pipeline de Normalización

### 1. Guardar conteos crudos (si no existen como layer)

```python
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()
```

### 2. Normalización por célula

```python
# Normalizar a 10,000 counts por célula (estándar)
sc.pp.normalize_total(adata, target_sum=1e4)
```

### 3. Log-transform

```python
sc.pp.log1p(adata)
```

### 4. Selección de HVG

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,           # estándar para datos humanos/ratón
    batch_key=None,             # si hay batch, pasarlo para regresión
    flavor="seurat_v3" if "counts" in adata.layers else "seurat",
    subset=False,               # no subsetear, solo marcar en .var
)
```

**Nota**: `seurat_v3` es más preciso pero requiere la capa `counts`. Si no existe, usar `seurat`.

Mostrar:

```python
sc.pl.highly_variable_genes(adata)
```

### 5. Opcional: Regresión de efectos de batch (si hay batch)

```python
if "batch" in adata.obs and adata.obs["batch"].nunique() > 1:
    # Comentarle al usuario que hay batch effect
    # Preguntar si quiere corregir: regress_out o combat
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mito"])
```

### 6. Escalado

```python
sc.pp.scale(adata, max_value=10)
```

## Parámetros Recomendados por Contexto

| Contexto | n_top_genes | target_sum | Notas |
|----------|-------------|------------|-------|
| 10x PBMC (humano) | 2000 | 1e4 | Estándar |
| Datos de ratón | 2000 | 1e4 | Similar a humano |
| Smart-seq2 | Todos | 1e6 | Full-length, no 10x |
| ATAC | 0 (no aplica) | N/A | Usar LSI (muon.atac) |
| Proteína (CITE) | 0 (no aplica) | N/A | CLR transform |

## Validación

Mostrar después de normalizar:

```python
# Verificar distribución
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts"])
print(f"Genes totales: {adata.n_vars}")
print(f"HVG seleccionados: {adata.var['highly_variable'].sum()}")
```

## Output

```python
{
    "stage": "normalize",
    "status": "completed",
    "summary": f"Normalizado. {n_hvg} HVG seleccionados de {n_vars} genes totales",
    "params": {
        "target_sum": 1e4,
        "log1p": True,
        "n_top_genes": 2000,
        "hvg_flavor": "seurat_v3",
        "scaled": True,
        "max_value": 10,
    },
    "n_hvg": int(adata.var["highly_variable"].sum()),
    "data_path": "/path/to/normalized_data.h5ad",
    "plots": ["hvg_plot.png"],
    "recommendations": [
        f"{n_hvg} HVG seleccionados → listo para clustering",
        "Sugiero correr PCA con n_comps=50",
    ]
}
```

## Hard Rules

- **NEVER delete, overwrite, or modify the original input file**. All normalized data goes to a NEW output path.
- Always save `counts` to `adata.layers` before normalizing
- If the user has protein data (CITE-Seq), do NOT use normalize_total — use muon's CLR instead
- If the user has ATAC data, redirect to `muon.atac.pp` (LSI instead of PCA)
- Do not subset the object to only HVG unless the user explicitly asks
- Save a copy before scaling (scaling is irreversible)
- The only person who removes original data is the human from the CLI.
