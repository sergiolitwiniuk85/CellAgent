# Proposal: Pipeline Multimodal Fixes

## Intent

Three methodological issues degrade CellAgent's multimodal pipeline: (1) unscaled feature concatenation causes PCA collapse, (2) overly aggressive HVG selection discards targeted panel biology, (3) missing spatial autocorrelation metrics.

## Scope

Pipeline stage: multimodal / spatial.

### In Scope
- Feature standardization before image-expression concatenation via `StandardScaler`
- Adaptive gene selection: all genes for targeted panels (<1000), top 2000 HVGs for whole transcriptome (>=1000)
- Spatial autocorrelation: Moran's I + Geary's C via `squidpy.gr.spatial_autocorr`

### Out of Scope
- GLCM optimization
- Compartment segmentation refactor
- Critic Agent implementation
- RAG-based annotation
- All Phase 2+ items

## Capabilities

### New Capabilities
None — these are methodological fixes to existing pipeline behavior.

### Modified Capabilities
None — implementation details change, not spec-level requirements.

## Approach

1. **StandardScaler**: Import `sklearn.preprocessing.StandardScaler` in both integration test scripts; `.fit_transform` image and expression features independently before `np.hstack`
2. **Adaptive genes**: Add `n_genes = adata.n_vars` check in normalize-agent; if < 1000, skip HVG selection; if >= 1000, select top 2000
3. **Spatial autocorr**: Insert `sq.gr.spatial_autocorr(adata, mode="moran")` and `mode="geary"` in spatial-agent pipeline; results stored in `adata.uns`

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `test_data/test_integration_fullfov_breast.py` | Modified | Add StandardScaler before concat |
| `test_data/test_integration_breast.py` | Modified | Add StandardScaler before concat |
| `skills/normalize-agent/SKILL.md` | Modified | Adaptive HVG selection logic |
| `skills/spatial-agent/SKILL.md` | Modified | Add Moran's I + Geary's C |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Feature scaling shifts downstream clusters | Medium | Compare cluster assignments before/after; document in report |
| All-gene mode may degrade clustering | Low | Targeted panels already sparse — more genes = more signal |

## Rollback Plan

Revert the 4 affected files to current state. No schema or data model changes — pure pipeline logic, reversible with `git checkout`. Old behavior preserved in git history.

## Dependencies

- `sklearn` (already in `environment.yml`)
- `squidpy` (already in `environment.yml`)

## Success Criteria

- [ ] PCA PC1 variance drops below 50% in integration tests after StandardScaler
- [ ] Targeted panel datasets use all genes (no HVG filter applied)
- [ ] Whole transcriptome datasets use top 2000 HVGs
- [ ] Moran's I and Geary's C values appear in `adata.uns` after spatial pipeline
