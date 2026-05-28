# Pipeline State — SCAI Shared State Schema

## State Object

Pipeline state is passed between agents as a dictionary/JSON with this structure:

```python
PipelineState = {
    # --- Identification ---
    "session_id": str,          # Unique analysis session ID
    "project_name": str,        # Project name (e.g. "pbmc_10k")
    
    # --- Data ---
    "data_path": str,           # Path to original data file
    "format": str,              # "h5ad" | "h5mu" | "zarr" | "h5ad+images"
    "modality": str,            # "rna" | "atac" | "protein" | "multimodal" | "spatial"
    
    # --- Main object (serialized to path) ---
    "current_data_path": str,   # Path to current .h5ad/.h5mu/.zarr file
    "current_data_type": str,   # "AnnData" | "MuData" | "SpatialData"
    
    # --- Dataset metadata ---
    "n_obs": int,               # Number of cells/observations
    "n_vars": int,              # Number of genes/features
    "layers": list[str],        # Available layers (counts, normalized, etc.)
    "obs_columns": list[str],   # Columns in .obs
    "var_columns": list[str],   # Columns in .var
    
    # --- Pipeline history ---
    "history": [
        {
            "stage": str,        # "data" | "qc" | "normalize" | "cluster" | "integration" | "spatial"
            "timestamp": str,    # ISO datetime
            "status": str,       # "completed" | "skipped" | "failed"
            "summary": str,      # What was done in this stage
            "params": dict,      # Parameters used
            "output_path": str,  # Path to generated file (if applicable)
        }
    ],
    
    # --- Per-stage results ---
    "qc_metrics": {
        "total_counts": float,
        "n_genes_by_counts": float,
        "pct_mito": float,
        "doublet_score": float | None,
        "filter_thresholds": dict,  # Thresholds used for filtering
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
        "cluster_key": str,         # Name in adata.obs (e.g. "leiden")
        "has_umap": bool,
        "has_markers": bool,
    },
    
    "integration": {
        "method": str | None,       # "WNN" | "MOFA+" | "concat" | None
        "modalities": list[str],    # Integrated modalities
    },
    
    "image_features": {
        "n_features": int,
        "n_images_processed": int,
        "levels_completed": list[str],  # ["morphology", "microenvironment", "compartments"]
        "feature_names": list[str],
        "feature_columns": list[str],   # column names added to .obs
        "plots": list[str],             # paths to generated plots
    },
    
    # --- Outputs ---
    "plots": list[dict],            # Paths to generated plots
    "report_path": str | None,      # Path to final report
}
```

## Lifecycle

1. **Start**: The orchestrator creates the state with `data_path` and what the user requested
2. **Each stage**: The agent reads current state, performs its work, updates the state
3. **Persistence**: After each stage, the orchestrator saves the state to engram
4. **Resume**: If the user returns later, the orchestrator retrieves the state and continues

## Rules

- Agents NEVER modify state directly in engram. They return the updated state to the orchestrator.
- The orchestrator is the sole persister.
- Always include `current_data_path` — it's the file the next agent will read.
- Plots are saved to `output/plots/{session_id}/` by default.
