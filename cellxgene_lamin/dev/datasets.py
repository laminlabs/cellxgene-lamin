def anndata_human_immune_cells():
    """Anndata object with semi-curated metadata."""
    import lamindb as ln

    adata = ln.core.datasets.anndata_human_immune_cells()
    adata.obs["sex_ontology_term_id"] = "PATO:0000384"
    adata.obs["organism"] = "human"
    adata.obs["sex"] = "unknown"
    # create some typos in the metadata
    adata.obs["tissue"] = adata.obs["tissue"].cat.rename_categories({"lung": "lungg"})
    # new donor ids
    adata.obs["donor"] = adata.obs["donor"].astype(str) + "-1"
    # drop animal cell
    adata = adata[adata.obs["cell_type"] != "animal cell", :]
    # remove columns that are reserved in the cellxgene schema
    adata.var.drop(columns=["feature_reference", "feature_biotype"], inplace=True)
    adata.raw.var.drop(
        columns=["feature_name", "feature_reference", "feature_biotype"], inplace=True
    )

    return adata
