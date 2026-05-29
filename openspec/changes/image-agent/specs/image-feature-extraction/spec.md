# Spec: Image Feature Extraction

**Capability**: `image-feature-extraction`
**Change**: `image-agent`
**Status**: Draft

## 1. Image Loading & Validation

- MUST load all image layers from `sdata.images` (Xenium, Visium, MERFISH, H&E, DAPI).
- MUST validate each image exists and has shape `(C, Y, X)` or `(Y, X)` before processing.
- MUST fail with a descriptive error if a referenced image key is missing or empty.
- MUST infer image type from `sdata` metadata — nucleic acid (DAPI) vs brightfield (H&E) vs fluorescence — and select appropriate feature extraction paths.
- SHOULD log image dimensions, dtype, and channel count on load.

## 2. Level 1 — Cell Morphology Features

- MUST extract per-cell morphology via `skimage.measure.regionprops` for each segmented mask.
- MUST compute: area, perimeter, eccentricity, solidity, equivalent diameter, and extent.
- MUST compute GLCM texture features (`skimage.feature.graycomatrix`) — contrast, dissimilarity, homogeneity, energy, correlation — aggregated per cell region.
- MUST return a DataFrame with one row per cell and columns prefixed by `morph_` and `texture_`.
- SHOULD skip cells smaller than a configurable `min_cell_area` (default 16 px²).

## 3. Level 2 — Microenvironment Features

- MUST compute local cell density as cell count within a configurable radius (default 50 µm) using `sklearn.neighbors.KDTree` or `scipy.spatial.cKDTree`.
- MUST compute distance from each cell centroid to the nearest tissue edge (binary mask boundary).
- MUST compute distance to nearest N neighbor centroids (N configurable, default 5).
- MUST return a DataFrame with columns prefixed by `microenv_`.
- SHOULD support multiple radii for density estimation via a `density_radii` parameter.

## 4. Level 3 — Tissue Compartments

- SHOULD segment tissue into compartments: epithelium, stroma, border, necrosis.
- SHOULD use Otsu thresholding on nuclear intensity followed by morphological closing to derive compartments when no annotation exists.
- SHOULD compute per-cell: distance to nearest compartment boundary, compartment label, and spatial gradient of compartment intensity.
- MUST fall back gracefully with a warning if compartment detection is unstable (e.g., uniform tissue).
- MUST NOT require pre-existing tissue annotations — compartments MUST be derived from image alone.

## 5. Multi-Image Merge

- MUST merge features from all processed images (H&E, DAPI, FOVs) into a single output DataFrame.
- MUST use spatial coordinates (centroid x, y in physical units) as the merge key across images.
- MUST handle one-to-many mappings (one cell → multiple FOVs) by averaging features and logging the count.
- MUST validate that coordinate systems are compatible before merging; MUST fail if units or registration differ.
- MUST prefix feature columns by source image key to prevent name collisions.

## 6. Pipeline Integration

- MUST store the merged feature matrix in `sdata["table"].obsm["image_features"]`.
- MUST store the following key metrics in `sdata["table"].obs`: `n_image_features`, `mean_cell_area`, `mean_eccentricity`, `tissue_compartment` (if Level 3 ran).
- SHOULD cast to dense numpy array for `.obsm` (DataFrame loses column names but preserves row order).
- MUST preserve cell ordering: row i in `obsm["image_features"]` corresponds to `sdata["table"].obs_names[i]`.

## 7. Output Generation

- MUST write a `features.csv` to the run output directory containing the merged feature DataFrame indexed by cell barcode.
- SHOULD generate QC plots: feature distribution histograms, PCA variance explained of image features, spatial scatter of top 2 PCs overlaid on tissue image.
- MUST NOT modify original image files in `sdata.images`.
- MUST NOT modify cell segmentation masks.

## 8. Report Integration

- The report-agent MUST include a dedicated image-derived features section when `sdata["table"].obsm["image_features"]` is present.
- The report MUST display a summary table: N cells, N features extracted, images processed, compartments found.
- The report MUST include the top-5 most variable features (by coefficient of variation) with per-cell scatter plots over spatial coordinates.
- The report SHOULD include PCA variance explained bar plot and spatial overlay of PC1.

## 9. Trace Integration

- The trace-agent MUST log every image processed with its key, shape, dtype, and feature count.
- The trace MUST record all parameters used: density radius, min cell area, number of neighbors, compartment method.
- The trace MUST record paths to all generated plots and `features.csv`.
- The trace SHOULD record timing per image and per feature level.
