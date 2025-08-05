import lamindb as ln
import bionty as bt
from django.db.models import Q
from cellxgene_lamin.dev import get_datasets_from_cxg

ln.track()

census_version = "2025-01-30"
previous_release = "2024-07-01"

cxg_datasets = get_datasets_from_cxg()
for idx, ds in enumerate(cxg_datasets):
    if idx % 10 == 0:
        print(f"annotating dataset {idx} of {len(cxg_datasets)}")
    af = ln.Artifact.filter(
        Q(key__contains=ds["dataset_id"]) & Q(key__contains=census_version)  # type: ignore
    ).one_or_none()
    if af is None:
        continue
    else:
        organism_ontology_ids = [
            organism["ontology_term_id"]  # type: ignore
            for organism in ds["organism"]  # type: ignore
        ]
        organism_records = bt.Organism.filter(
            ontology_id__in=organism_ontology_ids
        ).list()
        first_organism = organism_records[0]
        if first_organism.name == "house mouse":
            first_organism.name = "mouse"
        try:
            schema = ln.examples.cellxgene.get_cxg_schema(
                schema_version="5.2.0",
                field_types="ontology_id",
                organism=first_organism.name,
            )
        # An organism that we don't support first class
        except IndexError:
            pass

        curator = ln.curators.AnnDataCurator(af, schema)
        try:
            curator.validate()
            curator.save_artifact()

        except ln.errors.ValidationError as e:
            error_msg = str(e)
            # if the `tissue_type` is cell, there can be cell types in the `tissue_ontology_term_id`
            # we skip them for now
            if "not validated in feature 'tissue_ontology_term_id'" in error_msg:
                continue
            # some datasets have ethnicities like HANCESTRO:0005 || HANCESTRO:0008
            # We do not annotate them yet but potentially later
            elif (
                "term not validated in feature 'self_reported_ethnicity_ontology_term_id' in slot 'obs'"
                in error_msg
            ):
                continue
            # some datasets have diseases like MONDO:0005148 || MONDO:0008170
            # We do not annotate them yet but potentially later
            elif (
                "term not validated in feature 'disease_ontology_term_id' in slot 'obs'"
                in error_msg
            ):
                continue
            # there are datasets with missing columns like tissue_type
            elif "not in dataframe" in error_msg:
                continue
            else:
                raise

ln.finish()
