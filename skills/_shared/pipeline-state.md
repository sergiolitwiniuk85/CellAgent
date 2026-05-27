# Pipeline State — SCAI Shared State Schema

## Objeto de Estado

El estado del pipeline se pasa entre agentes como un diccionario/JSON con esta estructura:

```python
PipelineState = {
    # --- Identificación ---
    "session_id": str,          # ID único de la sesión de análisis
    "project_name": str,        # Nombre del proyecto (ej: "pbmc_10k")
    
    # --- Datos ---
    "data_path": str,           # Path al archivo de datos original
    "format": str,              # "h5ad" | "h5mu" | "zarr" | "h5ad+images"
    "modality": str,            # "rna" | "atac" | "protein" | "multimodal" | "spatial"
    
    # --- Objeto principal (serializado a path) ---
    "current_data_path": str,   # Path al archivo .h5ad/.h5mu/.zarr actual
    "current_data_type": str,   # "AnnData" | "MuData" | "SpatialData"
    
    # --- Metadata del dataset ---
    "n_obs": int,               # Número de células/observaciones
    "n_vars": int,              # Número de genes/features
    "layers": list[str],        # Layers disponibles (counts, normalized, etc.)
    "obs_columns": list[str],   # Columnas en .obs
    "var_columns": list[str],   # Columnas en .var
    
    # --- Pipeline history ---
    "history": [
        {
            "stage": str,        # "data" | "qc" | "normalize" | "cluster" | "integration" | "spatial"
            "timestamp": str,    # ISO datetime
            "status": str,       # "completed" | "skipped" | "failed"
            "summary": str,      # Resumen de lo que se hizo
            "params": dict,      # Parámetros usados
            "output_path": str,  # Path al archivo generado (si aplica)
        }
    ],
    
    # --- Resultados por etapa ---
    "qc_metrics": {
        "total_counts": float,
        "n_genes_by_counts": float,
        "pct_mito": float,
        "doublet_score": float | None,
        "filter_thresholds": dict,  # thresholds usados para filtrar
    },
    
    "normalization": {
        "method": str,              # "scanpy" | "scran" | "sctransform"
        "target_sum": float | None,
        "n_hvg": int,
        "hvg_method": str,          # "seurat_v3" | "seurat" | "cell_ranger"
    },
    
    "clustering": {
        "n_pcs": int,
        "resolution": float,
        "n_clusters": int,
        "cluster_key": str,         # nombre en adata.obs (ej: "leiden")
        "has_umap": bool,
        "has_markers": bool,
    },
    
    "integration": {
        "method": str | None,       # "WNN" | "MOFA+" | "concat" | None
        "modalities": list[str],    # modadalidades integradas
    },
    
    # --- Outputs ---
    "plots": list[dict],            # paths a plots generados
    "report_path": str | None,      # path al reporte final
}
```

## Ciclo de Vida

1. **Inicio**: El orquestador crea el estado con `data_path` y lo que el usuario dijo
2. **Cada etapa**: El agente lee el estado actual, hace su trabajo, actualiza el estado
3. **Persistencia**: Después de cada etapa, el orquestador guarda el estado en engram
4. **Continuación**: Si el usuario vuelve después, el orquestador recupera el estado y continúa

## Reglas

- Los agentes NUNCA modifican el estado directamente en engram. Devuelven el estado actualizado al orquestador.
- El orquestador es el único que persiste.
- Siempre incluir `current_data_path` — es el archivo que el próximo agente va a leer.
- Los plots se guardan en `output/plots/{session_id}/` por defecto.
