# Design: image-agent

## Technical Approach

Standalone agent that extracts quantitative morphological features from spatial transcriptomics tissue images using **classical computer vision** (scikit-image). Sits between spatial-agent and cluster-agent in the pipeline: spatial-agent renders and explores, image-agent quantifies, cluster-agent consumes expression + image features.

Strict **read-only** policy — never modifies images, masks, or `sdata["table"]` in place. Outputs go to a `features.csv` and `.obsm["image_features"]`.

## Architecture Decisions

### Standalone agent vs extending spatial-agent

`spatial-agent` focuses on Squidpy spatial statistics, neighborhood enrichment, and visualization — it's an analysis tool. Image feature extraction is a fundamentally different concern: image I/O, GPU-free CV pipelines, coordinate-space merging. Mixing them violates separation of concerns and would make both harder to maintain. A dedicated agent keeps each skill focused and independently testable.

### skimage over deep learning

| Criterion | skimage | Deep Learning |
|-----------|---------|---------------|
| GPU required | No | Yes |
| Deterministic | Yes | No (seed-dependent) |
| Scientifically validated | Decades | Months |
| Installation | pip/pip | CUDA + 2GB+ models |

For morphology (area, eccentricity, solidity, GLCM texture), skimage's `regionprops` is the gold standard. No neural network needed to measure a cell's area. Deep learning adds non-determinism, GPU dependency, and scientific reproducibility risk with zero benefit for these features.

### One merged CSV vs per-image outputs

Cells are the atomic unit of analysis. A single cell can appear across multiple images (H&E + DAPI + FOVs). Per-image outputs would force the consumer (cluster-agent) to merge by spatial coordinates anyway. Pre-merging once in image-agent is simpler and avoids repeated work. Merge key: physical centroid coordinates (x, y). One-to-many → average features, log count.

### CPU-only with optional cellpose

GPU dependency is a real barrier for bench scientists. The core pipeline uses skimage exclusively (CPU, lightweight). Cellpose is an optional add-on for users who want GPU segmentation — never required. This keeps the default path accessible to anyone with a laptop.

## Data Flow

```
SpatialData (images + masks + table)
    │
    ▼
image-agent (loads sdata, validates images)
    │
    ├── Level 1: skimage.measure.regionprops
    │   → morph_area, morph_eccentricity, morph_solidity, morph_perimeter,
    │     morph_extent, morph_equiv_diameter
    │   → skimage.feature.graycomatrix (per-cell GLCM)
    │   → texture_contrast, texture_dissimilarity, texture_homogeneity,
    │     texture_energy, texture_correlation
    │
    ├── Level 2: scipy.spatial.cKDTree
    │   → microenv_cell_density (within radius)
    │   → microenv_edge_distance (to mask boundary)
    │   → microenv_nndist_1..N (nearest neighbor distances)
    │
    ├── Level 3: skimage.filters.threshold_otsu + morphology
    │   → compartment_label (epithelium/stroma/border/necrosis)
    │   → compartment_boundary_distance
    │   → compartment_intensity_gradient
    │
    │  Merge by centroid (x, y) across all images
    │  Prefix columns by source image key
    │
    ▼
features.csv  +  sdata["table"].obsm["image_features"]
                       sdata["table"].obs["n_image_features"]
                       sdata["table"].obs["mean_cell_area"]
                       sdata["table"].obs["mean_eccentricity"]
                       sdata["table"].obs["tissue_compartment"]
    │
    ▼
cluster-agent (uses .X expression + .obsm["image_features"])
```

## File Changes

### New files

| File | Description |
|------|-------------|
| `skills/image-agent/SKILL.md` | Full agent instructions: 3-level extraction, merge, output spec |

### Modified files

| File | Change |
|------|--------|
| `skills/_shared/pipeline-state.md` | Add `image_features` section to state schema |
| `agents/scai-orchestrator/INSTRUCTIONS.md` | Add image-agent to pipeline diagram + delegation table |
| `skills/report-agent/SKILL.md` | Add Section 8: Image-Derived Features (summary table, top-5 variable features, spatial PC scatter) |
| `skills/trace-agent/SKILL.md` | Add image processing log stage — images processed, params, timing, plot paths |
| `environment.yml` | Add `scikit-image >= 0.22`, `tifffile`, `cellpose` (optional) |

## Interfaces

### Pipeline State Extension

```python
"image_features": {
    "n_features": int,
    "n_images_processed": int,
    "levels_completed": list[str],  # ["morphology", "microenvironment", "compartments"]
    "feature_names": list[str],
    "feature_columns": list[str],   # column names in .obs
    "plots": list[str],            # paths to generated plots
}
```

### Agent Output

```python
{
    "stage": "image-agent",
    "status": "completed",
    "n_images_processed": int,
    "n_features": int,
    "levels_completed": list[str],
    "features_csv": "output/{session_id}/features.csv",
    "summary": f"Extracted {n_features} features from {n_images} images",
    "plots": ["features_pca.png", "feature_distributions.png", "spatial_pc1.png"],
}
```

## Testing Strategy

No formal test framework exists in this project. Recommended manual validation:

1. **Xenium test data**: Run on a known Xenium dataset with nuclei masks. Verify morph_area matches known cell sizes (~50–200 µm² for human cells). GLCM texture values should be in [0, 1) range.
2. **Visium test data**: Run on Visium (spot-level). Verify features.csv has one row per spot, coordinates match `sdata["table"].obsm["spatial"]`.
3. **Reproducibility**: Run twice on same data → identical features (skimage is deterministic).
4. **CPU-only test**: Run with `cellpose` not installed → Level 1–3 work, no errors.
5. **Read-only verification**: `md5sum` original images before and after → identical.

## Rollback

- Remove `skills/image-agent/`
- Revert `pipeline-state.md`, orchestrator instructions, report-agent, trace-agent
- Revert `environment.yml`
