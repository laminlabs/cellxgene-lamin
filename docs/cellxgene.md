---
execute_via: python
---

# CELLxGENE: scRNA-seq

[CZ CELLxGENE](https://cellxgene.cziscience.com/) hosts one of the largest standardized collections of scRNA-seq data - LaminDB provides a streamlined interface to query and load it.

You can use the CELLxGENE data in two ways:

1. Query collections of `AnnData` objects.
2. Query slices from a single concatenated dataset without downloading everything via [TileDB-SOMA](https://github.com/single-cell-data/TileDB-SOMA).

To build similar data assets in-house:

1. See the [transfer guide](inv:docs#transfer) to zero-copy data to your own LaminDB instance.
2. See the [scRNA guide](inv:docs#scrna) to create a growing, standardized & versioned scRNA-seq dataset collection.

```{dropdown} Show me a screenshot

<img src="https://lamin-site-assets.s3.us-east-1.amazonaws.com/.lamindb/YHMYgXCfJTJvKPBm0000.png" width="700px">
```

```python
# pip install lamindb
!lamin init --modules bionty --storage ./test-cellxgene
```

```python
import lamindb as ln
```

Create the central query object for the public [`laminlabs/cellxgene`](https://lamin.ai/laminlabs/cellxgene) instance:

```python
db = ln.DB("laminlabs/cellxgene")
```

## Query for individual datasets

Every individual dataset in CELLxGENE is an `.h5ad` file that is stored as an artifact in LaminDB. Here is an exemplary query:

```python
users = db.User.lookup()
cell_types = db.bionty.CellType.lookup()

db.Artifact.filter(
    suffix=".h5ad",
    description__contains="immune",
    size__gt=1e9,  # size > 1GB
    cell_types__name__in=["B cell", "T cell"],  # cell types measured in AnnData
    created_by=users.sunnyosun,  # created by a specific user
).order_by("created_at").to_dataframe(
    include=["cell_types__name", "created_by__handle"]  # join with additional info
).head()
```

```{dropdown} What happens under the hood?

As you saw from inspecting `Artifact`, `Artifact.cell_types` relates artifacts with `bionty.CellType`.

The expression `cell_types__name__in` performs the join of the underlying registries and matches `bionty.CellType.name` to `["B cell", "T cell"]`.

Similar for `created_by`, which relates artifacts with `User`.
```

### Slice an individual dataset

Let's look at a CELLxGENE artifact and show its metadata using `.describe()`.

```python
artifact = db.Artifact.get(description="Mature kidney dataset: immune")
artifact.describe()
```

:::{dropdown} More ways of accessing metadata

Access just features:

```
artifact.features
```

Or get values associated with features:

```
artifact.features.get_values()
```

:::

To query & load a slice of the array data, you have several options:

1. Cache the artifact on disk and return the path to the cached data like: `artifact.cache() -> Path`
2. Cache & load the entire artifact into memory via `artifact.load() -> AnnData`
3. Stream the array using a (cloud-backed) accessor `artifact.open() -> AnnDataAccessor`

All of these option run much faster in the AWS `us-west-2` data center.

Cache:

```python
cache_path = artifact.cache()
cache_path
```

Cache & load:

```python
adata = artifact.load()
adata
```

Now we have an `AnnData` object, which stores observation annotations matching our artifact-level query in the `.obs` slot, and we can re-use almost the same query on the array-level.

```python
tissues = db.bionty.Tissue.lookup()
suspension_types = db.ULabel.filter(type__name="SuspensionType").lookup()
experimental_factors = db.bionty.ExperimentalFactor.lookup()

adata_slice = adata[
    adata.obs.cell_type.isin(
        [cell_types.dendritic_cell.name, cell_types.neutrophil.name]
    )
    & (adata.obs.tissue == tissues.kidney.name)
    & (adata.obs.suspension_type == suspension_types.cell.name)
    & (adata.obs.assay == experimental_factors.ln_10x_3_v2.name)
]
adata_slice
```

Stream, slice and load the slice into memory:

```python
adata_backed = artifact.open()
```

We now have an `AnnDataAccessor` object, which behaves much like an `AnnData`, and slicing looks similar to the query above.

See the slicing operation:

```python
adata_backed_slice = adata_backed[
    adata_backed.obs.cell_type.isin(
        [cell_types.dendritic_cell.name, cell_types.neutrophil.name]
    )
    & (adata_backed.obs.tissue == tissues.kidney.name)
    & (adata_backed.obs.suspension_type == suspension_types.cell.name)
    & (adata_backed.obs.assay == experimental_factors.ln_10x_3_v2.name)
]

adata_backed_slice.to_memory()
```

```python
adata_backed.close()
```

## Query collections of datasets

Let's search collections from CELLxGENE within the 2025-01-30 release: https://lamin.ai/laminlabs/cellxgene/collections, and then pick a top hit:

It's a Science paper and we can find more information on it using the [DOI](https://doi.org/10.1016/j.xgen.2023.100298) or CELLxGENE [collection id](https://cellxgene.cziscience.com/collections/af893e86-8e9f-41f1-a474-ef05359b1fb7). There are multiple versions of this collection.

```python
collection = db.Collection.get("quQDnLsMLkP3JRsC8gp5")
collection
```

```python
collection.versions.to_dataframe()
```

The collection groups artifacts.

```python
collection.artifacts.to_dataframe()
```

Let's now look at the collection that corresponds to a `cellxgene-census` release of `.h5ad` artifacts.

```python
collection = db.Collection.get(key="cellxgene-census", version="2025-01-30")
collection
```

You can query across artifacts by arbitrary metadata combinations, for instance:

```python
organisms = db.bionty.Organism.lookup()
experimental_factors = db.bionty.ExperimentalFactor.lookup()
tissues = db.bionty.Tissue.lookup()
suspension_types = db.ULabel.filter(type__name="SuspensionType").lookup()
features = db.Feature.lookup(
    return_field="name"
)  # here we choose to return .name directly
assays = db.bionty.ExperimentalFactor.lookup(return_field="name")
```

```python
query = collection.artifacts.filter(
    organisms=organisms.human,
    cell_types__in=[cell_types.dendritic_cell, cell_types.neutrophil],
    tissues=tissues.kidney,
    ulabels=suspension_types.cell,
    experimental_factors=experimental_factors.ln_10x_3_v2,
)
query = query.order_by("size")
query.to_dataframe().head()
```

### Slice a concatenated array

Let us now use the concatenated version of the `Census` collection: a `tiledbsoma` array that concatenates all `AnnData` arrays present in the `collection` we just explored. Slicing `tiledbsoma` works similar to slicing `DataFrame` or `AnnData`.

```python
value_filter = (
    f'{features.tissue} == "{tissues.brain.name}" and {features.cell_type} in'
    f' ["{cell_types.microglial_cell.name}", "{cell_types.neuron.name}"] and'
    f' {features.suspension_type} == "{suspension_types.cell.name}" and {features.assay} =='
    f' "{assays.ln_10x_3_v3}"'
)
value_filter
```

Query for the `tiledbsoma` array store that contains all concatenated expression data. It's a new dataset produced by concatenating all `AnnData`-like artifacts in the Census collection.

```python
census_artifact = db.Artifact.get(key="cell-census/2025-01-30/soma")
```

Run the slicing operation.

```python
human = "homo_sapiens"  # subset to human data

# open the array store for queries
with census_artifact.open() as store:
    # read SOMADataFrame as a slice
    cell_metadata = store["census_data"][human].obs.read(value_filter=value_filter)
    # concatenate results to pyarrow.Table
    cell_metadata = cell_metadata.concat()
    # convert to pandas.DataFrame
    cell_metadata = cell_metadata.to_pandas()

cell_metadata.head()
```

Create an `AnnData` object.

```python
from tiledbsoma import AxisQuery

with census_artifact.open() as store:
    experiment = store["census_data"][human]
    adata = experiment.axis_query(
        "RNA", obs_query=AxisQuery(value_filter=value_filter)
    ).to_anndata(
        X_name="raw",
        column_names={
            "obs": [
                features.assay,
                features.cell_type,
                features.tissue,
                features.disease,
                features.suspension_type,
            ]
        },
    )

adata.var = adata.var.set_index("feature_id")
adata
```

## Train ML models

You can either directly train ML models on very large collections of `AnnData`-like artifacts or on a single concatenated `tiledbsoma`-like artifact. For pros & cons of these approaches, see [this blog post](https://lamin.ai/blog/arrayloader-benchmarks).

### On a collection of arrays

{meth}`~lamindb.Collection.mapped` caches `AnnData` objects on disk and creates a map-style dataset that performs a virtual join of the features of the underlying `AnnData` objects.

```
from torch.utils.data import DataLoader

census_collection = db.Collection.get(name="cellxgene-census", version="2025-01-30")

dataset = census_collection.mapped(obs_keys=[features.cell_type], join="outer")

dataloader = DataLoader(dataset, batch_size=128, shuffle=True)

for batch in dataloader:
    pass

dataset.close()
```

For more background, see {doc}`docs:scrna-mappedcollection`.

### On a concatenated array

You can create streaming `PyTorch` dataloaders from `tiledbsoma` stores using `cellxgene_census` package.

```
import cellxgene_census.experimental.ml as census_ml

store = census_artifact.open()

experiment = store["census_data"][human]
experiment_datapipe = census_ml.ExperimentDataPipe(
    experiment,
    measurement_name="RNA",
    X_name="raw",
    obs_query=AxisQuery(value_filter=value_filter),
    obs_column_names=[features.cell_type],
    batch_size=128,
    shuffle=True,
    soma_chunk_size=10000,
)
experiment_dataloader = census_ml.experiment_dataloader(experiment_datapipe)

for batch in experiment_dataloader:
    pass

store.close()
```

For more background see [this guide](https://chanzuckerberg.github.io/cellxgene-census/notebooks/experimental/pytorch.html).

```{toctree}
:maxdepth: 1
:hidden:

cellxgene-curate
```
