from typing import Any

import argparse
import lamindb as ln
import bionty as bt
from django.db.models import Q
from cellxgene_lamin.dev import get_datasets_from_cxg, get_collections_from_cxg


parser = argparse.ArgumentParser()
parser.add_argument("--new", required=True, help="New census version")
parser.add_argument("--previous", required=True, help="Previous census version")
parser.add_argument("--track", action="store_true", help="Whether to track this run")
parser.add_argument(
    "--space", type=str, default=None, help="Space to use for registration"
)
parser.add_argument(
    "--limit",
    type=int,
    default=None,
    help="Limit number of datasets to process. Skips Collection & soma registration.",
)
args = parser.parse_args()

NEW_CENSUS_VERSION = args.new
PREVIOUS_CENSUS_VERSION = args.previous
CENSUS_S3PATH = f"s3://cellxgene-data-public/cell-census/{NEW_CENSUS_VERSION}/h5ads"

if args.track:
    track_kwargs: dict[str, Any] = {
        "params": {
            "new_census_version": NEW_CENSUS_VERSION,
            "previous_census_version": PREVIOUS_CENSUS_VERSION,
        },
    }
    if args.space is not None:
        track_kwargs["space"] = args.space
    ln.track(**track_kwargs)

# ---------------------------------------------------------------------------
# 1. Register artifacts from S3 & link to previous release
# ---------------------------------------------------------------------------

cxg_datasets: list[dict[str, Any]] = get_datasets_from_cxg()  # type: ignore
print(f"Found {len(cxg_datasets)} datasets from CELLxGENE API")

# Build lookup: dataset_id -> cxg metadata
cxg_lookup: dict[str, dict[str, Any]] = {ds["dataset_id"]: ds for ds in cxg_datasets}

artifacts_previous = ln.Artifact.filter(version_tag=PREVIOUS_CENSUS_VERSION).all()
print(f"Found {artifacts_previous.count()} artifacts from previous release")

h5ad_paths = list(ln.UPath(CENSUS_S3PATH).glob("*.h5ad"))
if args.limit is not None:
    h5ad_paths = h5ad_paths[: args.limit]

for path in h5ad_paths:
    dataset_id = path.stem

    # Check for previous version to pass to constructor
    artifact_previous = artifacts_previous.filter(
        key__endswith=f"{dataset_id}.h5ad"
    ).one_or_none()

    kwargs: dict[str, Any] = {}
    if artifact_previous is not None:
        kwargs["revises"] = artifact_previous

    artifact = ln.Artifact(path, **kwargs)
    artifact.version_tag = NEW_CENSUS_VERSION

    # Annotate with CXG metadata
    if dataset_id in cxg_lookup:
        artifact.n_observations = cxg_lookup[dataset_id]["cell_count"]
        artifact.description = cxg_lookup[dataset_id]["title"]

    artifact.save()

artifacts = ln.Artifact.filter(key__contains=NEW_CENSUS_VERSION).all()
print(f"Registered {len(artifacts)} artifacts")

if args.limit is None:
    # -------------------------------------------------------------------
    # 2. Register collections
    # -------------------------------------------------------------------

    collection = ln.Collection(
        artifacts,
        key="cellxgene-census",
        revises=ln.Collection.filter(
            key="cellxgene-census", version_tag=PREVIOUS_CENSUS_VERSION
        ).one(),
    )
    collection.version_tag = NEW_CENSUS_VERSION
    collection.save()

    cxg_collections: list[dict[str, Any]] = get_collections_from_cxg()  # type: ignore
    ln.settings.creation.search_names = False

    for collection_meta in cxg_collections:
        keys = [
            f"cell-census/{NEW_CENSUS_VERSION}/h5ads/{dataset['dataset_id']}.h5ad"
            for dataset in collection_meta["datasets"]
        ]
        collection_artifacts = artifacts.filter(key__in=keys).all()
        if collection_artifacts.count() > 0:
            previous_collection = ln.Collection.filter(
                reference=collection_meta["collection_id"],
                version_tag=PREVIOUS_CENSUS_VERSION,
            ).one_or_none()

            collection_kwargs: dict[str, Any] = {
                "key": collection_meta["name"],
                "description": collection_meta["doi"],
                "reference": collection_meta["collection_id"],
                "reference_type": "CELLxGENE Collection ID",
            }
            if previous_collection is not None:
                collection_kwargs["revises"] = previous_collection

            collection_record = ln.Collection(collection_artifacts, **collection_kwargs)
            collection_record.version_tag = NEW_CENSUS_VERSION
            if collection_record._state.adding:
                collection_record.save()

    ln.settings.creation.search_names = True

    # -------------------------------------------------------------------
    # 3. Register the soma store
    # -------------------------------------------------------------------

    soma_path = f"s3://cellxgene-data-public/cell-census/{NEW_CENSUS_VERSION}/soma"
    previous_soma = ln.Artifact.filter(
        description=f"Census {PREVIOUS_CENSUS_VERSION}"
    ).one()
    soma_artifact = ln.Artifact(
        soma_path,
        description=f"Census {NEW_CENSUS_VERSION}",
        revises=previous_soma,
    )
    soma_artifact.version_tag = NEW_CENSUS_VERSION
    soma_artifact.save()
    print(soma_artifact)

# ---------------------------------------------------------------------------
# 4. Annotate artifacts (validate & curate)
# ---------------------------------------------------------------------------

cxg_datasets_to_annotate: list[dict[str, Any]] = cxg_datasets
if args.limit is not None:
    cxg_datasets_to_annotate = cxg_datasets_to_annotate[: args.limit]

for idx, ds in enumerate(cxg_datasets_to_annotate):
    if idx % 10 == 0:
        print(f"Annotating dataset {idx} of {len(cxg_datasets_to_annotate)}")

    af = ln.Artifact.filter(
        Q(key__contains=ds["dataset_id"]) & Q(key__contains=NEW_CENSUS_VERSION)
    ).one_or_none()
    if af is None:
        continue

    organism_ontology_ids = [
        organism["ontology_term_id"] for organism in ds["organism"]
    ]
    organism_records = bt.Organism.filter(
        ontology_id__in=organism_ontology_ids
    ).to_list()
    first_organism = organism_records[0]
    if first_organism.name == "house mouse":
        first_organism.name = "mouse"

    try:
        schema = ln.examples.cellxgene.create_cellxgene_schema(
            field_types="ontology_id",
            organism=first_organism.name,
        )
    except IndexError:
        continue

    curator = ln.curators.AnnDataCurator(af, schema)
    try:
        curator.validate()
        curator.save_artifact()
    except ln.errors.ValidationError as e:
        error_msg = str(e)
        if "not validated in feature 'tissue_ontology_term_id'" in error_msg:
            continue
        elif (
            "term not validated in feature 'self_reported_ethnicity_ontology_term_id' in slot 'obs'"
            in error_msg
        ):
            continue
        elif (
            "term not validated in feature 'disease_ontology_term_id' in slot 'obs'"
            in error_msg
        ):
            continue
        elif "not in dataframe" in error_msg:
            continue
        else:
            raise

if args.track:
    ln.finish()
