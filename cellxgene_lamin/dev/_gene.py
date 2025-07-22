import pandas as pd
from rich import print


def register_genes() -> None:
    """Note that this bulk saves genes which can lead to IntegrityErrors."""
    import bionty as bt
    import lamindb as ln

    bt.Organism.from_source(
        name="synthetic construct", source=bt.Source.get("4tsksCMX")
    ).save()
    bt.Organism.from_source(name="sars-2", source=bt.Source.get("4tsksCMX")).save()

    organisms = bt.Organism.lookup(field=bt.Organism.scientific_name)
    genes_files = {
        "homo_sapiens": "https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/genes_homo_sapiens.csv.gz",
        "mus_musculus": "https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/genes_mus_musculus.csv.gz",
        "synthetic_construct": "https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/genes_ercc.csv.gz",
        "severe_acute_respiratory_syndrome_coronavirus_2": "https://github.com/chanzuckerberg/single-cell-curation/raw/main/cellxgene_schema_cli/cellxgene_schema/gencode_files/genes_sars_cov_2.csv.gz",
    }

    # register all genes for each organism
    for organism_name, genes_file in genes_files.items():
        print(f"[bold orange]registering {organism_name} genes")
        organism_record = getattr(organisms, organism_name)
        if organism_name == "synthetic_construct":
            df = pd.read_csv(
                genes_file,
                header=None,
                index_col=0,
                names=["ensembl_gene_id", "symbol", "chromosome", "length", "biotype"],
            )
            gene_records = bt.Gene.from_values(
                df.index,
                field=bt.Gene.ensembl_gene_id,
                organism=organism_record,
                create=True,
            )
            for gene in gene_records:
                gene.symbol = df.loc[gene.ensembl_gene_id].symbol
                gene.save()
            ln.save([r for r in gene_records if r._state.adding])
        elif organism_name == "severe_acute_respiratory_syndrome_coronavirus_2":
            # currently broken https://github.com/chanzuckerberg/single-cell-curation/issues/1415
            pass
        else:
            df = pd.read_csv(genes_file, header=None, index_col=0)
            gene_records = bt.Gene.from_values(
                df.index, field=bt.Gene.ensembl_gene_id, organism=organism_record
            )
            ln.save([r for r in gene_records if r._state.adding])
        validated = bt.Gene.validate(
            df.index, field=bt.Gene.ensembl_gene_id, organism=organism_record
        )
        # register legacy genes manually
        if organism_name == "severe_acute_respiratory_syndrome_coronavirus_2":
            pass
        else:
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
