def anndata_human_immune_cells(populate_registries):
    """Anndata object with semi-curated metadata."""
    import lamindb as ln

    adata = ln.dev.datasets.anndata_human_immune_cells()
    adata.obs["sex_ontology_term_id"] = "PATO:0000384"
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

    if populate_registries:
        import bionty as bt

        verbosity = ln.settings.verbosity
        ln.settings.verbosity = "error"
        cell_types = bt.CellType.from_values(adata.obs["cell_type"])
        ln.save(cell_types)
        organism = bt.Organism.from_public(name="human")
        organism.save()
        genes = bt.Gene.from_values(
            adata.var_names, field=bt.Gene.ensembl_gene_id, organism=organism
        )
        ln.save(genes[:-10])
        ln.settings.verbosity = verbosity

    return adata
