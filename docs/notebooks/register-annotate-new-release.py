from typing import Any
import argparse
import lamindb as ln
import bionty as bt
from lamin_utils import logger
from django.db.models import Q
from django.core.exceptions import ObjectDoesNotExist
from cellxgene_lamin.dev._cxg_rest_api import (
    get_datasets_from_cxg,
    get_collections_from_cxg,
)

ln.settings.sync_git_repo = "https://github.com/laminlabs/cellxgene-lamin"


parser = argparse.ArgumentParser()
parser.add_argument("--new", required=True, help="New census version")
parser.add_argument(
    "--previous",
    required=False,
    default=None,
    help="Previous census version (full release only)",
)
parser.add_argument("--track", action="store_true", help="Whether to track this run")
parser.add_argument(
    "--space", type=str, default=None, help="Space to use for registration"
)
parser.add_argument(
    "--pre-release",
    action="store_true",
    dest="pre_release",
    help="Register datasets not yet in LaminDB as pre-release (labeled, not version-tagged).",
)
parser.add_argument(
    "--smoke",
    action="store_true",
    help="Limits number of datasets to process to 2. "
    "Skips Collection & soma registration.",
)

args = parser.parse_args()
logger.info(
    f"starting run | new={args.new} | previous={args.previous} | smoke={args.smoke} "
    f"| track={args.track} | pre_release={args.pre_release}"
)

NEW_CENSUS_VERSION = args.new
PREVIOUS_CENSUS_VERSION = args.previous
CENSUS_S3PATH = f"s3://cellxgene-data-public/cell-census/{NEW_CENSUS_VERSION}/h5ads"

if args.track:
    track_kwargs: dict[str, Any] = {
        "params": {
            "new_census_version": NEW_CENSUS_VERSION,
            "previous_census_version": PREVIOUS_CENSUS_VERSION,
            "pre_release": args.pre_release,
        }
    }
    if args.space is not None:
        track_kwargs["space"] = args.space
    ln.track("Rrq1bb328HH4", **track_kwargs)
    logger.info(f"tracking enabled | space={args.space}")


if args.smoke:
    ln.examples.cellxgene.save_cellxgene_defaults()
    logger.info("smoke mode: saved cellxgene defaults")


cxg_datasets: list[dict[str, Any]] = get_datasets_from_cxg()  # type: ignore
logger.info(f"found {len(cxg_datasets)} datasets from CellxGene")

cxg_lookup: dict[str, dict[str, Any]] = {ds["dataset_id"]: ds for ds in cxg_datasets}  # type: ignore

registered_ids: set[str] = set()
registered_artifacts: dict[str, Any] = {}

if args.pre_release:
    import cellxgene_census as cxc

    # -----------------------------------------------------------------------
    # 1. Register pre-release artifacts (datasets not already in LaminDB / LTS)
    # -----------------------------------------------------------------------

    # get or create the pre-release label
    pre_release_label = ln.ULabel.filter(name="pre-release").one_or_none()
    if pre_release_label is None:
        pre_release_label = ln.ULabel(name="pre-release").save()
    logger.info("pre-release ULabel ready")

    # cleanup: delete broken pre-release artifacts from previous runs
    existing_pre_release = ln.Artifact.filter(ulabels=pre_release_label)
    logger.info(
        f"checking {existing_pre_release.count()} existing pre-release artifacts for breakage"
    )
    for af in existing_pre_release:
        try:
            af.open()
        except Exception as e:
            logger.warning(f"broken, deleting | key={af.key} | {e}")
            af.delete(permanent=True)

    # register new datasets
    for dataset_id, ds in cxg_lookup.items():
        # skip if already registered in LaminDB (already in LTS or a prior run)
        if (
            ln.Artifact.filter(key__endswith=f"{dataset_id}.h5ad").one_or_none()
            is not None
        ):
            logger.info(f"already registered, skipping | dataset_id={dataset_id}")
            continue

        # new datasets are NOT in the versioned S3 path — resolve via latest census
        try:
            uri = cxc.get_source_h5ad_uri(dataset_id, census_version="latest")["uri"]
        except Exception as e:
            logger.warning(
                f"could not resolve h5ad uri | dataset_id={dataset_id} | {e}"
            )
            continue

        artifact = ln.Artifact(uri, description=ds["title"])

        # check the file is accessible before registering
        try:
            artifact.open()
        except Exception as e:
            logger.warning(f"broken dataset, skipping | dataset_id={dataset_id} | {e}")
            continue

        artifact.n_observations = ds["cell_count"]
        artifact.save()
        artifact.ulabels.add(pre_release_label)

        registered_ids.add(dataset_id)
        registered_artifacts[dataset_id] = artifact
        logger.info(f"registered pre-release | dataset_id={dataset_id}")

        if args.smoke and len(registered_ids) >= 2:
            break

    logger.info(f"registered {len(registered_ids)} pre-release artifacts")

    # -----------------------------------------------------------------------
    # 2. Register pre-release collections (only collections NOT already in LTS)
    # -----------------------------------------------------------------------
    if registered_artifacts and not args.smoke:
        cxg_collections: list[dict[str, Any]] = get_collections_from_cxg()  # type: ignore
        logger.info(f"found {len(cxg_collections)} CellxGene collections")

        ln.settings.creation.search_names = False
        for collection_meta in cxg_collections:
            # skip if this collection already exists in LTS (a version-tagged record)
            in_lts = ln.Collection.filter(
                reference=collection_meta["collection_id"],
                version_tag__isnull=False,
            ).exists()
            if in_lts:
                logger.info(
                    f"collection already in LTS, skipping | id={collection_meta['collection_id']}"
                )
                continue

            # brand-new collection → all its registered datasets are pre-release
            member_artifacts = [
                registered_artifacts[d["dataset_id"]]
                for d in collection_meta["datasets"]
                if d["dataset_id"] in registered_artifacts
            ]
            if not member_artifacts:
                continue

            # chain weekly pre-release reruns of the same collection
            previous_collection = ln.Collection.filter(
                reference=collection_meta["collection_id"],
                ulabels=pre_release_label,
            ).one_or_none()

            collection_kwargs: dict[str, Any] = {
                "key": f"pre-release/{collection_meta['name']}",
                "description": collection_meta["doi"],
                "reference": collection_meta["collection_id"],
                "reference_type": "CELLxGENE Collection ID",
            }
            if previous_collection is not None:
                collection_kwargs["revises"] = previous_collection
                logger.info(
                    f"revising pre-release collection: {collection_meta['name']}"
                )
            else:
                logger.info(
                    f"creating pre-release collection: {collection_meta['name']}"
                )

            collection_record = ln.Collection(member_artifacts, **collection_kwargs)
            collection_record.save()
            collection_record.ulabels.add(pre_release_label)

        ln.settings.creation.search_names = True

else:
    # -----------------------------------------------------------------------
    # 1. Register artifacts (full LTS release)
    # -----------------------------------------------------------------------
    h5ad_paths = list(ln.UPath(CENSUS_S3PATH).glob("*.h5ad"))
    logger.info(f"found {len(h5ad_paths)} h5ad paths in {CENSUS_S3PATH}")
    if args.smoke:
        h5ad_paths = h5ad_paths[:2]
        logger.info("smoke mode: limiting to 2 h5ad paths")

    for h5ad_path in h5ad_paths:
        dataset_id = h5ad_path.stem
        registered_ids.add(dataset_id)
        artifact_previous = (
            ln.Artifact.filter(
                key__endswith=f"{dataset_id}.h5ad",
                key__contains=PREVIOUS_CENSUS_VERSION,
            ).one_or_none()
            if PREVIOUS_CENSUS_VERSION
            else None
        )
        kwargs: dict[str, Any] = {}
        if artifact_previous is not None:
            kwargs["revises"] = artifact_previous
            logger.info(f"revising existing artifact for dataset_id={dataset_id}")
        else:
            logger.info(f"registering new artifact for dataset_id={dataset_id}")

        artifact = ln.Artifact(h5ad_path, **kwargs)
        artifact.version_tag = NEW_CENSUS_VERSION
        if dataset_id in cxg_lookup:
            artifact.description = cxg_lookup[dataset_id]["title"]
            artifact.n_observations = cxg_lookup[dataset_id]["cell_count"]
        artifact.save()
        registered_artifacts[dataset_id] = artifact

    new_afs = ln.Artifact.filter(key__contains=NEW_CENSUS_VERSION)
    logger.info(
        f"registered {len(registered_ids)} artifacts for census version {NEW_CENSUS_VERSION}"
    )

    if not args.smoke:
        # -------------------------------------------------------------------
        # 2. Register collections
        # -------------------------------------------------------------------
        logger.info("registering top-level cellxgene-census collection")
        collection = ln.Collection(
            new_afs,
            key="cellxgene-census",
            revises=ln.Collection.filter(
                key="cellxgene-census", version_tag=PREVIOUS_CENSUS_VERSION
            ).one(),
        )
        collection.version_tag = NEW_CENSUS_VERSION
        collection.save()
        logger.info("saved top-level collection")

        cxg_collections = get_collections_from_cxg()  # type: ignore
        logger.info(f"found {len(cxg_collections)} CellxGene collections")

        ln.settings.creation.search_names = False
        for collection_meta in cxg_collections:
            keys = [
                f"cell-census/{NEW_CENSUS_VERSION}/h5ads/{dataset['dataset_id']}.h5ad"
                for dataset in collection_meta["datasets"]
            ]
            collection_artifacts = new_afs.filter(key__in=keys)
            if collection_artifacts.count() > 0:
                previous_collection = ln.Collection.filter(
                    reference=collection_meta["collection_id"],
                    version_tag=PREVIOUS_CENSUS_VERSION,
                ).one_or_none()

                collection_kwargs = {
                    "key": collection_meta["name"],
                    "description": collection_meta["doi"],
                    "reference": collection_meta["collection_id"],
                    "reference_type": "CELLxGENE Collection ID",
                }
                if previous_collection is not None:
                    collection_kwargs["revises"] = previous_collection
                    logger.info(
                        f"revising collection: {collection_meta['name']} (id={collection_meta['collection_id']})"
                    )
                else:
                    logger.info(
                        f"creating new collection: {collection_meta['name']} (id={collection_meta['collection_id']})"
                    )

                collection_record = ln.Collection(
                    collection_artifacts, **collection_kwargs
                )
                collection_record.version_tag = NEW_CENSUS_VERSION
                collection_record.save()
            else:
                logger.warning(
                    f"no matching artifacts for collection: {collection_meta['name']} (id={collection_meta['collection_id']}), skipping"
                )

        ln.settings.creation.search_names = True

        # -------------------------------------------------------------------
        # 3. Register the soma store
        # -------------------------------------------------------------------
        logger.info("registering soma store")
        soma_path = f"s3://cellxgene-data-public/cell-census/{NEW_CENSUS_VERSION}/soma"
        previous_soma = ln.Artifact.filter(
            description=f"Census {PREVIOUS_CENSUS_VERSION}"
        ).one()
        new_soma_af = ln.Artifact(
            soma_path,
            description=f"Census {NEW_CENSUS_VERSION}",
            revises=previous_soma,
        )
        new_soma_af.version_tag = NEW_CENSUS_VERSION
        new_soma_af.save()
        logger.info(f"saved soma artifact: {new_soma_af}")

# ---------------------------------------------------------------------------
# 4. Annotate artifacts (validate & curate)
# ---------------------------------------------------------------------------
logger.info("starting annotation of artifacts")

for idx, ds in enumerate(cxg_datasets):
    if ds["dataset_id"] not in registered_ids:
        continue
    if idx % 10 == 0:
        logger.info(f"annotating dataset {idx} of {len(cxg_datasets)}")

    if args.pre_release:
        af = ln.Artifact.filter(
            ulabels=pre_release_label,
            key__endswith=f"{ds['dataset_id']}.h5ad",
        ).one_or_none()
    else:
        af = ln.Artifact.filter(
            Q(key__contains=ds["dataset_id"]) & Q(key__contains=NEW_CENSUS_VERSION)
        ).one_or_none()

    if af is None:
        logger.warning(f"no artifact found for dataset_id={ds['dataset_id']}, skipping")
        continue

    organism_ontology_ids = [
        organism["ontology_term_id"] for organism in ds["organism"]
    ]
    organism_records = bt.Organism.filter(
        ontology_id__in=organism_ontology_ids
    ).to_list()

    if not organism_records:
        logger.warning(
            f"skipping dataset_id={ds['dataset_id']}: no organism records found"
        )
        continue
    first_organism = organism_records[0]
    if first_organism.name == "house mouse":
        first_organism.name = "mouse"

    try:
        schema = ln.examples.cellxgene.create_cellxgene_schema(
            field_types="ontology_id", organism=first_organism.name
        )
    except ObjectDoesNotExist:
        logger.warning(
            f"skipping dataset_id={ds['dataset_id']}: bt.Source not found for "
            f"organism={first_organism.name}. "
            f"run bt.Gene.add_source(organism='{first_organism.name}') to fix."
        )
        continue
    except IndexError:
        logger.warning(
            f"skipping dataset_id={ds['dataset_id']}: IndexError while creating "
            f"schema for organism={first_organism.name}"
        )
        continue

    curator = ln.curators.AnnDataCurator(af, schema)

    try:
        curator.validate()
        curator.save_artifact()
        logger.info(f"successfully validated and saved dataset_id={ds['dataset_id']}")

    except ln.errors.ValidationError as e:
        error_msg = str(e)
        if "not validated in feature 'tissue_ontology_term_id'" in error_msg:
            logger.warning(
                f"skipping dataset_id={ds['dataset_id']}: tissue_ontology_term_id not validated"
            )
            continue
        elif (
            "term not validated in feature 'self_reported_ethnicity_ontology_term_id' in slot 'obs'"
            in error_msg
        ):
            logger.warning(
                f"skipping dataset_id={ds['dataset_id']}: self_reported_ethnicity_ontology_term_id not validated"
            )
            continue
        elif (
            "not validated in feature 'disease_ontology_term_id' in slot 'obs'"
            in error_msg
        ):
            logger.warning(
                f"skipping dataset_id={ds['dataset_id']}: disease_ontology_term_id not validated"
            )
            continue
        elif "not in dataframe" in error_msg:
            logger.warning(
                f"skipping dataset_id={ds['dataset_id']}: feature not in dataframe"
            )
            continue
        elif "no Organism found in source for the given" in error_msg:
            logger.warning(
                f"skipping dataset_id={ds['dataset_id']}: ontology not in dataframe"
            )
            continue
        else:
            logger.error(
                f"unhandled ValidationError for dataset_id={ds['dataset_id']}: {error_msg}"
            )
            continue

if args.track:
    ln.finish()
