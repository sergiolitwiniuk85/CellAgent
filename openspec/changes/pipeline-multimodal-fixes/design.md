# Design: Pipeline Multimodal Fixes

## Technical Approach

Three independent methodological fixes with zero coupling. Each fix is a localized change to one pipeline phase — no refactoring, no new abstractions.

| Fix | Phase | Change | Insertion Point |
|-----|-------|--------|----------------|
| Feature Standardization | 4 — Shared kNN | `StandardScaler` per modality before `np.hstack` | After smoothed matrices, before `np.hstack` |
| Adaptive Gene Selection | 2 — Alignment | Conditional HVG count based on `n_genes` | Replace hardcoded `[-20:]` |
| Spatial Autocorrelation | 7 — Leiden Clustering | `sq.gr.spatial_autocorr` Moran's I + Geary's C | After `sc.tl.leiden`, before summary |

## Architecture Decisions

### Fix 1: Feature Standardization

| Decision | Choice | Alternatives | Rationale |
|----------|--------|-------------|-----------|
| Scaler type | `StandardScaler` (z-score) | `MinMaxScaler`, `RobustScaler` | PCA assumes unit variance; z-score is standard preprocessing for PCA. `MinMaxScaler` doesn't handle outliers. `RobustScaler` is overkill — image features have bounded ranges, expression is sparse but not heavy-tailed. |
| Scaling scope | Per modality separately | Global scaling, per-feature scaling | Global scaling would let the larger modality (image: ~25 features vs expr: 20-2000) still dominate. Per-modality ensures equal contribution. |
| Transform timing | After kNN smoothing | Before kNN, after hstack | Smoothing preserves scale; scaling after smoothing is correct — we standardize the final fused representation per modality. |
| State persistence | Save scaler objects | Discard after fit | Required for reproducibility and potential inverse_transform in reporting. Saved alongside output. |

### Fix 2: Adaptive Gene Selection

| Decision | Choice | Alternatives | Rationale |
|----------|--------|-------------|-----------|
| Threshold | `n_genes < 1000` | `n_genes < 500`, `n_genes < 2000` | Xenium panels are 250–500 genes, Visium/CosMx WTA are 5000–20000. 1000 cleanly separates targeted panels from whole transcriptome. |
| Panel behavior | Use ALL genes | Use top N (e.g., 100) | If only 300 genes are measured, selecting 20 discards 93% of the data. All-gene mode preserves full panel signal. |
| WTA behavior | Top 2000 HVGs | Top 5000, top 1000 | 2000 is standard for scRNA-seq HVG selection (Seurat default). Matches normalize-agent's `n_top_genes=2000`. |

### Fix 3: Spatial Autocorrelation

| Decision | Choice | Alternatives | Rationale |
|----------|--------|-------------|-----------|
| Method | `sq.gr.spatial_autocorr` | Manual implementation, `esda` | Squidpy is already a dependency; wraps `esda` internally. Single function call for both Moran's I and Geary's C. |
| Metrics | Both Moran's I + Geary's C | Only Moran's I | Moran's I detects global autocorrelation; Geary's C is more sensitive to local variation. Complementary — report both. |
| Scope | Full FOV only (has Leiden) | Also add Leiden to crop version | Crop test is a lighter integration test. Adding Leiden + autocorr would change its purpose. Keep crop focused on core pipeline validation. |

## Data Flow

### Phase 2 (modified for Fix 2)
```
expr_crop (n_cells × n_genes)
  → if n_genes < 1000: top_var_genes_idx = all genes
  → else: gene_var = np.var(expr_crop, axis=0); top_var_genes_idx = top 2000
```

### Phase 4 (modified for Fix 1)
```
image_feat_smoothed (n_cells × n_img_feats)
  → StandardScaler().fit_transform() → img_scaled
expr_smoothed (n_cells × n_expr_feats)
  → StandardScaler().fit_transform() → expr_scaled
  → np.hstack([img_scaled, expr_scaled]) → multimodal_matrix
  → PCA(...)
```

### Phase 7 (modified for Fix 3)
```
adata with leiden clusters
  → sq.gr.spatial_autocorr(adata, mode="moran")
  → sq.gr.spatial_autocorr(adata, mode="geary")
  → adata.uns["moranI"], adata.uns["gearyC"] populated
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `test_data/test_integration_fullfov_breast.py` | Modify | Fix 1: import StandardScaler, scale before hstack, save scalers. Fix 2: adaptive gene count. Fix 3: add Moran's I + Geary's C after Leiden. Update summary.json to include n_expression_genes. |
| `test_data/test_integration_breast.py` | Modify | Fix 1: same scaling logic. Fix 2: same adaptive gene logic. No Fix 3 (no Leiden phase in crop version). |
| `skills/normalize-agent/SKILL.md` | Modify | Add HVG Rule: if total genes < 1000, skip HVG selection (targeted panel). Document 1000-threshold in Recommended Parameters table. |
| `skills/spatial-agent/SKILL.md` | Modify | Add Section 6: Spatial Autocorrelation with Moran's I + Geary's C after clustering. Include code example and output structure. |

## Testing Strategy

| Fix | Verification |
|-----|-------------|
| Fix 1 | Run both scripts. PCA PC1 variance should drop vs baseline (previously inflated by unscaled morph_area dominance). `scaler_img` and `scaler_expr` files saved to output. |
| Fix 2 | Run with Xenium data (n_genes < 1000): expression features in output should equal total gene count. Run with Visium data (n_genes >= 1000): should report 2000 HVGs. Log line prints adaptive count. |
| Fix 3 | Run fullfov script. Verify `adata.uns['moranI']` and `adata.uns['gearyC']` are non-empty DataFrames with expected columns. Summary includes spatial autocorrelation metrics. |

## Open Questions

- Should scaler parameters be saved for reproducibility? **Yes**. Save both scaler objects via `pickle` or `joblib` alongside scaler position in data flow. On re-run, skip fit and use saved scaler for inverse transform in reporting.
- Should Moran's I threshold for "spatially variable" be configurable? **Not yet**. Current design stores raw statistics without thresholding. Configurable p-value threshold (`moran_p_threshold`) can be added in Phase 2 if reporting needs "top N spatially variable genes" filtering.
