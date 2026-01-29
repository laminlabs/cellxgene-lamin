---
execute_via: python
---

<!-- #region -->

# Curate `AnnData` based on the CELLxGENE schema

This guide shows how to curate an AnnData object against the [CELLxGENE schema v5.2.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.2.0/schema.md).

````{admonition} Summary
To ingest validate & annotated datasets adhering to a CELLxGENE Schema, call

```bash
!cellxgene-schema --version # should print 5.2.0
!cellxgene-schema validate small_cxg_curated.h5ad  # validation
```

using a shell, and then

```python
schema = ln.examples.cellxgene.create_cellxgene_schema(version="5.2.0")
ln.Artifact("…", schema=schema).save()  # annotation (re-validates ontologies, but not some other details)
```
````

<!-- #endregion -->

```python
# pip install lamindb pronto
# cellxgene-schema has pinned dependencies. Therefore we recommend installing it into a separate environment using `uv` or `pipx`
# uv tool install cellxgene-schema==5.2.3

!lamin init --storage ./test-cellxgene-curate --modules bionty
```

```python
import lamindb as ln
import bionty as bt

ln.track()
```

## The CELLxGENE schema

As a first step, we generate the specific CELLxGENE schema which adds missing sources to the instance:

```python
cxg_schema = ln.examples.cellxgene.create_cellxgene_schema("5.2.0")
```

```python
cxg_schema.describe()
```

The schema has two components:

```python
cxg_schema.slots["var"].describe()
```

```python
cxg_schema.slots["obs"].describe()
```

In the following, we will validate a dataset the CELLxGENE schema and curate it.

## Validate and curate metadata

Let's start with an AnnData object that we would like to curate.
We are writing it to disk to run [CZI's cellxgene-schema CLI tool](https://github.com/chanzuckerberg/single-cell-curation) which verifies whether an on-disk h5ad dataset adheres all requirements of CELLxGENE including the CELLxGENE schema.

```python
adata = ln.examples.datasets.small_dataset3_cellxgene(
    with_obs_typo=True, with_var_typo=True
)
adata.write_h5ad("small_cxg.h5ad")
adata
```

Initially, the `cellxgene-schema` validator of CZI does not pass and we need to curate the dataset.

```python
!MPLBACKEND=agg uvx cellxgene-schema validate small_cxg.h5ad
```

CELLxGENE requires all observations to be annotated.
If information for a specific column like `disease_ontology_term_id` is not available, CELLxGENE requires to fall back to default values like "normal" or "unknown".
Let's save these defaults to the instance using {func}`lamindb.examples.cellxgene.save_cellxgene_defaults`:

```python
ln.examples.cellxgene.save_cellxgene_defaults()
```

Now we can start curating the dataset:

```python
curator = ln.curators.AnnDataCurator(adata, cxg_schema)
try:
    curator.validate()
except ln.errors.ValidationError:
    pass
```

The error shows invalid genes are present in the dataset.
Let's remove them from both the `adata` and `adata.raw` objects:

```python
adata = adata[
    :, ~adata.var.index.isin(curator.slots["var"].cat.non_validated["index"])
].copy()
if adata.raw is not None:
    raw_data = adata.raw.to_adata()
    raw_data = raw_data[
        :, ~raw_data.var.index.isin(curator.slots["var"].cat.non_validated["index"])
    ].copy()
    adata.raw = raw_data
```

As we've subsetted the AnnData object, we have to recreate the `AnnDataCurator` to validate again:

```python
curator = ln.curators.AnnDataCurator(adata, cxg_schema)
try:
    curator.validate()
except ln.errors.ValidationError as e:
    print(e)
```

The validation error tells us that we're missing several columns.
The reason is simple:
CELLxGENE requires all `obs` metadata to be stored as ontology IDs in `entity_ontology_term_id` columns.
Therefore, we first translate the `name` based `obs` columns into the required format.

```python
adata.obs
```

```python
# Add missing assay column
adata.obs["assay_ontology_term_id"] = "EFO:0005684"
# Add `entity_ontology_term_id` columns by translating names to ontology IDs
standardization_map = {
    "self_reported_ethnicity": (
        bt.Ethnicity,
        "self_reported_ethnicity_ontology_term_id",
    ),
    "cell_type": (bt.CellType, "cell_type_ontology_term_id"),
}

for col, (bt_class, new_col) in standardization_map.items():
    adata.obs[new_col] = bt_class.standardize(
        adata.obs[col], field="name", return_field="ontology_id"
    )
# Drop the name columns because CELLxGENE disallows them
adata.obs = adata.obs.drop(columns=list(standardization_map.keys()))
```

```python
try:
    curator.validate()
except ln.errors.ValidationError:
    pass
```

An error is shown for the tissue label “UBERON:0002048XXX” because it contains a few extra `X` - a typo.
Let’s fix it:

```python
adata.obs["tissue_ontology_term_id"] = adata.obs[
    "tissue_ontology_term_id"
].cat.rename_categories({"UBERON:0002048XXX": "UBERON:0002048"})
```

Now `validate` should pass.

```python
# recreate the AnnDataCurator to refresh cached categoricals
curator = ln.curators.AnnDataCurator(adata, cxg_schema)
curator.validate()
```

## Save artifact

We can now save the curated artifact:

```python
artifact = curator.save_artifact(key="examples/dataset-curated-against-cxg.h5ad")
```

```python
artifact.describe()
```

## Validating using cellxgene-schema

To validate the now curated AnnData object using [CZI's cellxgene-schema CLI tool](https://github.com/chanzuckerberg/single-cell-curation), we need to write the AnnData object to disk.

```python
adata.write("small_cxg_curated.h5ad")
```

```python
# %%bash -e
!MPLBACKEND=agg uvx cellxgene-schema validate small_cxg_curated.h5ad
```

```{note}

The CELLxGENE Schema is designed to validate all metadata for adherence to ontologies.
It does not reimplement all rules of the cellxgene schema and we therefore recommend running the [cellxgene-schema](https://github.com/chanzuckerberg/single-cell-curation) if full adherence beyond metadata is a necessity.
```
