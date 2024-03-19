from typing import TYPE_CHECKING, Iterable

from lamindb import Artifact
from lnschema_core import Registry


def register_organisms(cxg_datasets):
    import bionty as bt
    import lamindb as ln

    # register all organisms
    ncbitaxon_source = bt.BiontySource.filter(source="ncbitaxon").one()

    organisms_meta = set()
    for cxg_dataset in cxg_datasets:
        organisms_meta.update({i["ontology_term_id"] for i in cxg_dataset["organism"]})

    organisms_records = bt.Organism.from_values(
        organisms_meta, field=bt.Organism.ontology_id, bionty_source=ncbitaxon_source
    )
    # rename 'house mouse' to 'mouse'
    for r in organisms_records:
        if r.name == "house mouse":
            r.name = "mouse"
    ln.save(organisms_records, parents=False)


def annotate_organisms(artifacts: Artifact, cxg_datasets: Iterable):
    import bionty as bt
    import lamindb as ln

    feature_organism = ln.Feature.filter(name="organism").one()

    for cxg_dataset in cxg_datasets:
        # get registered file record based on dataset_id
        artifact = artifacts.filter(
            key__contains=cxg_dataset["dataset_id"]
        ).one_or_none()
        if artifact is None:
            continue

        # annotate artifacts with organisms
        organism_ontology_ids = [i["ontology_term_id"] for i in cxg_dataset["organism"]]
        organism_records = bt.Organism.filter(
            ontology_id__in=organism_ontology_ids
        ).list()
        artifact.labels.add(organism_records, feature=feature_organism)
