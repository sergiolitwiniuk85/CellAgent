# Tasks: image-agent

## Review Workload Forecast

| Field | Value |
|-------|-------|
| Estimated changed lines | ~280 (200 new + 80 modified) |
| 400-line budget risk | Low |
| Chained PRs recommended | No |
| Delivery strategy | ask-always |
| Chain strategy | size-exception |

Decision needed before apply: Yes
Chained PRs recommended: No
Chain strategy: size-exception
400-line budget risk: Low

**Note**: 280 lines is well under the 400-line budget. Single PR. `size-exception` because the new SKILL.md (~150 lines) is an LLM-instruction artifact that reads denser than code. No chaining needed.

## Phase 1: Foundation

- [x] 1.1 Create `skills/image-agent/SKILL.md` — full agent instructions covering Level 1 (regionprops + GLCM), Level 2 (cKDTree density + edge + NND), Level 3 (Otsu compartments + gradients), multi-image merge by centroid (x,y), output to `features.csv` + `.obsm["image_features"]` + `.obs` columns, QC plot generation, CPU-only default with optional cellpose
- [x] 1.2 Update `skills/_shared/pipeline-state.md` — add `image_features` block to state schema: `n_features`, `n_images_processed`, `levels_completed`, `feature_names`, `feature_columns`, `plots`
- [x] 1.3 Update `environment.yml` — add `scikit-image >= 0.22`, `tifffile`, `cellpose` (commented as optional)

## Phase 2: Pipeline Wiring

- [x] 2.1 Update `agents/scai-orchestrator/INSTRUCTIONS.md` pipeline diagram — add `image-agent` arrow between spatial-agent and cluster-agent
- [x] 2.2 Update `agents/scai-orchestrator/INSTRUCTIONS.md` delegation table — add row for image-agent: `skills/image-agent/SKILL.md` | `general`

## Phase 3: Report & Trace Updates

- [x] 3.1 Update `skills/report-agent/SKILL.md` — add Section: Image-Derived Features (summary table, top-5 most variable features with spatial scatter, PCA variance bar plot, PC1 spatial overlay)
- [x] 3.2 Update `skills/trace-agent/SKILL.md` — add image processing log stage (images processed, params used, plot paths, timing per image per level)

## Phase 4: Verification

- [x] 4.1 Manual validation — run end-to-end pipeline on Xenium test data; verify `features.csv` row count matches cells, `obsm["image_features"]` shape, `.obs` columns populated, and original images unmodified
