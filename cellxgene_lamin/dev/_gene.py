import pandas as pd
from lamin_utils import logger


def register_genes() -> None:
    """Note that this bulk saves genes which can lead to IntegrityErrors."""
    import bionty as bt
    import lamindb as ln

    genes_files = {
        "homo_sapiens": "https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/genes_homo_sapiens.csv.gz",
        "mus_musculus": "https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/genes_mus_musculus.csv.gz",
        "synthetic_construct": "https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/genes_ercc.csv.gz",
        "severe_acute_respiratory_syndrome_coronavirus_2": "https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/genes_sars_cov_2.csv.gz",
    }

    # Register missing Organism & Source if necessary
    bt.Organism.from_source(
        name="synthetic construct", source=bt.Source.get("4tsksCMX")
    ).save()
    bt.Organism.from_source(name="sars-2", source=bt.Source.get("4tsksCMX")).save()

    synthetic_df = pd.read_csv(
        genes_files["synthetic_construct"],
        header=None,
        index_col=0,
        names=["ensembl_gene_id", "name", "chromosome", "length", "biotype"],
    )
    source = bt.Source(
        entity="bionty.Gene",
        name="gencode_ercc",
        version="1.0.0",
        organism="synthetic construct",
        currently_used=False,
        description="GENCODE with ERCC spike-ins",
        url=genes_files["synthetic_construct"],
    ).save()
    bt.Gene.add_source(source=source, df=synthetic_df)

    sars_2_df = pd.read_csv(
        genes_files["severe_acute_respiratory_syndrome_coronavirus_2"],
        header=None,
        index_col=0,
        names=["ensembl_gene_id", "name", "chromosome", "length", "biotype"],
    )
    source = bt.Source(
        entity="bionty.Gene",
        name="sars_cov_2",
        version="1.0.0",
        organism="sars-2",
        currently_used=False,
        description="GENCODE of Sars Cov 2",
        url=genes_files.get("severe_acute_respiratory_syndrome_coronavirus_2"),
    ).save()
    bt.Gene.add_source(source=source, df=sars_2_df)

    organisms = bt.Organism.lookup(field=bt.Organism.scientific_name)

    # register all genes for each organism
    for organism_name, genes_file in genes_files.items():
        logger.info(f"CELLXGENE: registering {organism_name} genes")
        organism_record = getattr(organisms, organism_name)
        if organism_record.name == "house mouse":
            organism_record.name = "mouse"

        df = pd.read_csv(genes_file, header=None, index_col=0)
        gene_records = bt.Gene.from_values(
            df.index, field=bt.Gene.ensembl_gene_id, organism=organism_record
        )
        ln.save([r for r in gene_records if r._state.adding])
        validated = bt.Gene.validate(
            df.index, field=bt.Gene.ensembl_gene_id, organism=organism_record
        )
        # register legacy genes manually
        new_records = []
        for gene_id in df.index[~validated]:
            new_records.append(
                bt.Gene(
                    ensembl_gene_id=gene_id,
                    symbol=df.loc[gene_id][1],
                    organism=organism_record,
                )
            )
        ln.save(new_records)
