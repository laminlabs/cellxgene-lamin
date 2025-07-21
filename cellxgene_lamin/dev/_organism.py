from collections.abc import Iterable

from lamindb import Artifact


def register_organisms(cxg_datasets: Iterable) -> None:
    import bionty as bt
    import lamindb as ln

    # register all organisms
    ncbitaxon_source = bt.Source.filter(name="ncbitaxon").first()

    organisms_ontology_ids = set()
    for cxg_dataset in cxg_datasets:
        organisms_ontology_ids.update(
            {organism["ontology_term_id"] for organism in cxg_dataset["organism"]}
        )

    organisms_records = bt.Organism.from_values(
        organisms_ontology_ids, field=bt.Organism.ontology_id, source=ncbitaxon_source
    )
    # rename 'house mouse' to 'mouse'
    for record in organisms_records:
        if record.name == "house mouse":
            record.name = "mouse"
    ln.save(organisms_records)


def curate_organisms(artifacts: Artifact, cxg_datasets: Iterable) -> None:
    import bionty as bt
    import lamindb as ln

    for cxg_dataset in cxg_datasets:
        artifact = artifacts.filter(
            key__contains=cxg_dataset["dataset_id"]
        ).one_or_none()
        if artifact is None:
            continue

        organism_ontology_ids = [
            organism["ontology_term_id"] for organism in cxg_dataset["organism"]
        ]
        organism_records = bt.Organism.filter(
            ontology_id__in=organism_ontology_ids
        ).list()
        artifact.organisms.set(organism_records)
