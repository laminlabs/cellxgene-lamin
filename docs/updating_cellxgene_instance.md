# CELLxGENE Dataset Ingestion Runbook

## Purpose

Update the `laminlabs/cellxgene` instance whenever CZ CELLxGENE publishes a new LTS Census data release.

## Trigger

A new entry appears on the [LTS Census data releases page](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html#list-of-lts-census-data-releases).

**Current latest LTS**: `2025-11-08` (Census schema 2.4.0, Dataset schema 7.0.0, 1845 datasets)

## Steps

### 1. Check schema changes

1. Fetch the new release's **Dataset schema version** from the LTS release page.
2. Read the full schema spec at `https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/<VERSION>/schema.md`.
3. Diff against the previous schema version. Look for:
   - New or removed `obs`/`var`/`uns` fields.
   - New ontology sources or changed ontology versions.
   - Changed cardinality or delimiter rules (e.g. multi-value `disease` with `' || '` in 7.0.0).
   - New organism support.
   - New embedding types.
4. If changes are found, update `create_cellxgene_schema()` and/or `save_cellxgene_defaults()` in [`lamindb/examples/cellxgene/_cellxgene.py`](https://github.com/laminlabs/lamindb/blob/main/lamindb/examples/cellxgene/_cellxgene.py):
   - Add/remove entries in `categoricals_to_spec`.
   - Update `CategorySpec` field references, defaults, or `needs_organism` flags.
   - Update `save_cellxgene_defaults()` if new default/control values are required.
   - Update the `CELLxGENEOrganisms` literal if new organisms are supported.
5. Open a PR against `laminlabs/lamindb`, get review, merge.

### 2. Update instance versions

1. Open the transform notebook at [lamin.ai/laminlabs/cellxgene/transform/SsAKLtg3y83R0000](https://lamin.ai/laminlabs/cellxgene/transform/SsAKLtg3y83R0000).
2. Update version strings / parameters to match the new LTS release.
3. Re-run the notebook end-to-end.
4. Verify the instance metadata reflects the new schema version.

### 3. Registration & annotation script

The script lives at [github.com/laminlabs/cellxgene](cellxgene-lamin/docs/notebooks/register-annotate-new-release.py).
It handles: registering all datasets as Artifacts with revision chains, registering Collections, creating the census soma store, and annotating with schema.
Commit to cellxgene-lamin repo.

### 4. Test registration (smoke test)

Run with --limit 1 against a local LaminDB instance:

```bash
python register-annotate-new-release.py --new 2025-11-08 --previous 2025-01-30 --limit 1
```

Verify: schema validation passes, all `_ontology_term_id` features are populated, the artifact otype is AnnData.

5. Run full registration on AWS SageMaker

Launch a SageMaker instance with sufficient resources.
Execute the full run:

```bash
python register-annotate-new-release.py --new 2025-11-08 --previous 2025-01-30 --track --space cellxgene
```

Expected runtime: ~24 hours+.
Monitor for failures; log any datasets that error out for manual follow-up.

### 6. Post-run verification

Run the following checks against the `laminlabs/cellxgene` instance:

| Check                                                      | Expected                                                                     |
| ---------------------------------------------------------- | ---------------------------------------------------------------------------- |
| Total artifacts with `otype="AnnData"`                     | Matches number of datasets in the new LTS release (e.g. 1845 for 2025-11-08) |
| All AnnData artifacts have `_ontology_term_id` annotations | 100% coverage                                                                |
| All AnnData artifacts have references attached             | 100% coverage                                                                |
| New SOMA store exists for this release                     | Present and accessible                                                       |
| No duplicate registrations from prior releases             | Confirmed                                                                    |

```python
import lamindb as ln
import bionty as bt
import pandas as pd

EXPECTED_DATASET_COUNT = 1845  # update per LTS release

# --- 1. Total count ---
artifacts = ln.Artifact.filter(otype="AnnData")
n_anndata = artifacts.count()
assert n_anndata == EXPECTED_DATASET_COUNT, f"Expected {EXPECTED_DATASET_COUNT}, got {n_anndata}"

# --- 2. Ontology ID coverage (every artifact must have all required annotations) ---
required_links = [
    "cell_types",       # CellType
    "diseases",         # Disease
    "tissues",          # Tissue
    "experimental_factors",  # ExperimentalFactor (assay)
    "developmental_stages",
    "phenotypes",       # sex
    "ethnicities",
    "organisms",
]

missing: dict[str, list[str]] = {}
for af in artifacts.all():
    for link in required_links:
        if not getattr(af, link, None) or getattr(af, link).count() == 0:
            missing.setdefault(link, []).append(af.uid)

if missing:
    summary = {k: len(v) for k, v in missing.items()}
    raise AssertionError(f"Artifacts missing ontology links:\n{summary}\nFirst UIDs per link: { {k: v[:3] for k, v in missing.items()} }")

# --- 3. References ---
no_ref = [af.uid for af in artifacts.all() if af.references.count() == 0]
assert len(no_ref) == 0, f"{len(no_ref)} artifacts have no references. First 5: {no_ref[:5]}"

# --- 4. Schema compliance (spot-check a sample) ---
no_schema = [af.uid for af in artifacts.all() if af.schema is None]
assert len(no_schema) == 0, f"{len(no_schema)} artifacts have no schema. First 5: {no_schema[:5]}"

# --- 5. Version lineage (new artifacts should be new versions of previously registered ones) ---
new_artifacts = artifacts.order_by("-created_at")[:EXPECTED_DATASET_COUNT]
no_lineage = [
    af.uid for af in new_artifacts
    if af.version is not None and af.stem_uid is not None and af.previous_versions.count() == 0
]
# Artifacts that are v1 (no previous versions expected) are fine; flag those claiming a version but having no predecessors
assert len(no_lineage) == 0, f"{len(no_lineage)} versioned artifacts have no previous versions. First 5: {no_lineage[:5]}"

# --- 6. No duplicates ---
keys = list(artifacts.values_list("key", flat=True))
dupes = [k for k in set(keys) if keys.count(k) > 1]
assert len(dupes) == 0, f"Duplicate artifact keys found: {dupes[:10]}"

# --- 7. SOMA store exists ---
soma_artifacts = ln.Artifact.filter(otype="tiledbsoma").order_by("-created_at")
assert soma_artifacts.count() > 0, "No SOMA store found"
latest_soma = soma_artifacts.first()
print(f"Latest SOMA store: {latest_soma.uid} created {latest_soma.created_at}")
```

## Rollback

If verification fails:

1. Identify the failing datasets and check logs.
2. If schema generation is the issue, fix `_cellxgene.py` and re-run from Step 2.
3. If a small number of datasets fail, re-run registration for just those datasets.

## Key URLs

- **LTS releases**: https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html#list-of-lts-census-data-releases
- **Dataset schema specs**: https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema
- **Schema generation code**: https://github.com/laminlabs/lamindb/blob/main/lamindb/examples/cellxgene/_cellxgene.py
- **Instance version notebook**: https://lamin.ai/laminlabs/cellxgene/transform/SsAKLtg3y83R0000
- **Registration notebooks**: https://github.com/laminlabs/cellxgene-lamin/tree/main/docs/notebooks
