# Proposal: image-agent

## Intent

Spatial transcriptomics datasets (Xenium, Visium, MERFISH) include tissue images rich in morphology. The spatial-agent renders them but never extracts features. We need an agent that quantifies cell morphology, microenvironments, and compartment gradients from images and feeds them into clustering.

## Scope

### In Scope
- `skills/image-agent/SKILL.md` with full agent instructions
- Level 1: Cell morphology — area, eccentricity, solidity, GLCM texture via skimage
- Level 2: Microenvironment — local density, edge distance, nearest neighbor
- Level 3: Tissue compartments + gradients (epithelium/stroma/border/necrosis)
- `features.csv` merged by spatial coords across all images
- Attach to `sdata["table"].obsm["image_features"]` + key metrics in `.obs`
- Read-only — strict "never modify" policy
- Pipeline: spatial → **image-agent** → cluster
- Update **report-agent** with image-derived features section (summary table, top features, spatial scatter plots)
- Update **trace-agent** with image processing log (images processed, features extracted, parameters used, plots)
- Pipeline state schema extension (`image_features` section)

### Out of Scope
- Image registration (SpatialData handles it)
- Cell segmentation (uses existing masks)
- Pathology-grade classification
- GPU deep learning (cellpose is optional)

## Capabilities

### New Capabilities
- `image-feature-extraction`: Cell morphology, microenvironment, and tissue compartment features from spatial tissue images

### Modified Capabilities
- None

## Approach

1. Create `skills/image-agent/SKILL.md` — loads SpatialData, extracts 3 levels
2. Merge multi-image features by spatial coordinates
3. Write `features.csv`, attach to `.obsm` + `.obs`
4. Update `pipeline-state.md`, orchestrator instructions, `environment.yml`

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `skills/image-agent/SKILL.md` | New | Full agent instructions |
| `skills/_shared/pipeline-state.md` | Modified | Add `image_features` section |
| `agents/scai-orchestrator/INSTRUCTIONS.md` | Modified | Pipeline diagram + delegation |
| `skills/report-agent/SKILL.md` | Modified | Add image-derived features section |
| `skills/trace-agent/SKILL.md` | Modified | Add image processing log |
| `environment.yml` | Modified | scikit-image, tifffile, cellpose |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| GPU dependency (cellpose) | Medium | CPU-only fallback via skimage |
| Xenium scale (100k+ cells) | Medium | Batch per FOV, memory-mapped reads |
| Visium registration | Low | Handled by SpatialData |

## Rollback Plan

- Remove `skills/image-agent/` directory
- Revert `pipeline-state.md` and orchestrator instructions
- Revert `environment.yml` dependency additions
- Pipeline returns to original data → qc → normalize → spatial → cluster → report → trace order

## Dependencies

- `scikit-image >= 0.22` — morphology, GLCM texture, segmentation utilities
- `tifffile` — multi-page TIFF support
- `cellpose` (optional) — GPU-accelerated segmentation fallback
- Existing: SpatialData, scanpy, numpy, pandas

## Success Criteria

- [ ] `skills/image-agent/SKILL.md` exists and passes markdown lint
- [ ] Agent extracts Level 1–3 features from Xenium test data (nuclei masks + H&E)
- [ ] `features.csv` merges multiple images (H&E + DAPI) by spatial coordinates
- [ ] Features land in `sdata["table"].obsm["image_features"]` + `.obs`
- [ ] Pipeline runs end-to-end: data → qc → normalize → spatial → image-agent → cluster
- [ ] CPU-only path works without cellpose installed
- [ ] Original image files are unmodified after run
