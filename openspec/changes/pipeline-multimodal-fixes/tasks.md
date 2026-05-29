# Tasks: Pipeline Multimodal Fixes

## Review Workload Forecast

| Fix | Files | Estimated Δ Lines |
|-----|-------|------------------|
| Fix 1 — StandardScaler | `test_integration_fullfov_breast.py`, `test_integration_breast.py` | ~15 per script = 30 |
| Fix 2 — Adaptive gene selection | Both scripts + `normalize-agent/SKILL.md` | ~20 per script + 10 SKILL.md = ~50 |
| Fix 3 — Spatial autocorrelation | `test_integration_fullfov_breast.py` + `spatial-agent/SKILL.md` | ~10 + ~15 = ~25 |
| **Total** | **4 files** | **~105 lines** |

Well under 400 lines → single PR. No pre-apply decisions needed.

---

## Phase 1: Feature Standardization (Fix 1)

- [x] **1.1** — Import `StandardScaler` from `sklearn.preprocessing` in both scripts
- [x] **1.2** — In `test_integration_fullfov_breast.py` Phase 4: `.fit_transform` `image_feat_smoothed` and `expr_smoothed` separately before `np.hstack`
- [x] **1.3** — Same change in `test_integration_breast.py` Phase 4 (identical insertion point)
- [x] **1.4** — Save both scaler objects via `joblib.dump` to `OUT / "scaler_img.pkl"` and `OUT / "scaler_expr.pkl"`
- [x] **1.5** — Update summary.json to include scaler paths and log a note about per-modality scaling

**Check**: PCA PC1 variance drops vs baseline (no longer dominated by unscaled `morph_area`).

---

## Phase 2: Adaptive Gene Selection (Fix 2)

- [x] **2.1** — In both scripts' Phase 2: read `n_genes` from expression matrix, wrap HVG selection in `if n_genes >= 1000: top 2000` / `else: use all genes`
- [x] **2.2** — Replace hardcoded `[-20:]` top-var slice with computed `n_top_var = min(2000 if n_genes >= 1000 else n_genes, 2000)`
- [x] **2.3** — Update summary print and `summary.json` `n_expression_genes` to reflect dynamic count
- [x] **2.4** — In `normalize-agent/SKILL.md`: add row to "Recommended Parameters by Context" table for targeted panels (`n_genes < 1000`: all genes, skip HVG). Add Hard Rule: "If `adata.n_vars < 1000`: use all genes, do NOT call `highly_variable_genes`"

**Check**: Xenium run uses all genes (300-500). With synthetic whole-transcriptome data (≥1000), selects 2000 HVGs.

---

## Phase 3: Spatial Autocorrelation (Fix 3)

- [x] **3.1** — In `test_integration_fullfov_breast.py` Phase 7: after `sc.tl.leiden`, add `sq.gr.spatial_autocorr(adata, mode="moran")` and `mode="geary"` (spatial coordinates already in `adata.obsm['spatial']`)
- [x] **3.2** — Add print logging: show top-5 Moran's I and Geary's C genes
- [x] **3.3** — Append spatial autocorrelation keys to `summary.json`
- [x] **3.4** — In `spatial-agent/SKILL.md`: add Section 6 "Spatial Autocorrelation" after clustering, with code example (`sq.gr.spatial_autocorr`), output structure (`adata.uns['moranI']`, `.uns['gearyC']`), and interpretation notes
- [x] **3.5** — No changes to `test_integration_breast.py` (no Leiden phase)

**Check**: `adata.uns['moranI']` and `adata.uns['gearyC']` are non-empty DataFrames with `I`, `pval_sim`, `pval_z_sim`, `var_n`, `var_z` columns.

---

## Phase 4: Verification

- [x] **4.1** — Run `test_integration_breast.py`: verify no errors, check PCA variance drop, confirm `n_expression_genes` matches total genes (targeted panel)
- [x] **4.2** — Run `test_integration_fullfov_breast.py`: same checks, plus verify `adata.uns['moranI']` / `['gearyC']` populated
- [x] **4.3** — If CSV golden files exist, update them; if not, verify output tables are well-formed
- [x] **4.4** — `git diff --stat` to confirm ~105 lines total
