from __future__ import annotations

from typing import Any

import bionty as bt
import click
import lamindb as ln
from django.core.exceptions import ObjectDoesNotExist
from django.db.models import Q
from lamin_utils import logger

from cellxgene_lamin.dev._cxg_rest_api import (
    get_collections_from_cxg,
    get_datasets_from_cxg,
)

ln.settings.sync_git_repo = "https://github.com/laminlabs/cellxgene-lamin"


# ---------------------------------------------------------------------------
# Shared helpers (also tested directly)
# ---------------------------------------------------------------------------


def cleanup_broken_pre_release_artifacts(pre_release_label: Any) -> None:
    """Delete any existing pre-release artifacts that are no longer accessible."""
    existing = ln.Artifact.filter(ulabels=pre_release_label)
    logger.info(
        f"checking {existing.count()} existing pre-release artifacts for breakage"
    )
    for af in existing:
        try:
            af.open()
        except Exception as e:
            logger.warning(f"broken, deleting | key={af.key} | {e}")
            af.delete(permanent=True)


def register_pre_release_artifacts(
    cxg_lookup: dict[str, dict[str, Any]],
    pre_release_label: Any,
    cxc: Any,
    smoke: bool = False,
) -> tuple[set[str], dict[str, Any]]:
    """Register datasets from the latest census not already in LaminDB."""
    registered_ids: set[str] = set()
    registered_artifacts: dict[str, Any] = {}

    existing_dataset_ids = {
        key.rsplit("/", 1)[-1].removesuffix(".h5ad")
        for key in ln.Artifact.filter(key__endswith=".h5ad").values_list(
            "key", flat=True
        )
    }
    logger.info(f"found {len(existing_dataset_ids)} already-registered datasets")

    # one bulk census read: get all census dataset IDs and their h5ad URIs
    with cxc.open_soma(census_version="latest") as census:
        h5ad_base = census.uri.replace("/soma", "/h5ads").rstrip("/")
        census_df = (
            census["census_info"]["datasets"]
            .read()
            .concat()
            .to_pandas()[["dataset_id", "dataset_h5ad_path"]]
        )
    census_uri_map: dict[str, str] = {
        row["dataset_id"]: f"{h5ad_base}/{row['dataset_h5ad_path']}"
        for _, row in census_df.iterrows()
    }
    logger.info(f"found {len(census_uri_map)} datasets in census")

    cxg_lookup_filtered = {
        dataset_id: ds
        for dataset_id, ds in cxg_lookup.items()
        if dataset_id not in existing_dataset_ids and dataset_id in census_uri_map
    }
    logger.info(f"{len(cxg_lookup_filtered)} datasets to register")

    items = list(cxg_lookup_filtered.items())
    if smoke:
        items = items[:2]

    for dataset_id, ds in items:
        uri = census_uri_map[dataset_id]
        artifact = ln.Artifact(uri, description=ds["title"])
        artifact.n_observations = ds["cell_count"]
        artifact.save()
        artifact.ulabels.add(pre_release_label)

        registered_ids.add(dataset_id)
        registered_artifacts[dataset_id] = artifact
        logger.info(f"registered pre-release | dataset_id={dataset_id}")

    return registered_ids, registered_artifacts


def register_pre_release_collections(
    cxg_collections: list[dict[str, Any]],
    registered_artifacts: dict[str, Any],
    pre_release_label: Any,
) -> None:
    """Register CellxGene collections that are not yet in any LTS release."""
    lts_collection_ids = set(
        ln.Collection.filter(version_tag__isnull=False).values_list(
            "reference", flat=True
        )
    )
    existing_pre_release = {
        c.reference: c for c in ln.Collection.filter(ulabels=pre_release_label)
    }

    ln.settings.creation.search_names = False
    for collection_meta in cxg_collections:
        if collection_meta["collection_id"] in lts_collection_ids:
            logger.info(
                f"collection already in LTS, skipping | id={collection_meta['collection_id']}"
            )
            continue

        member_artifacts = [
            registered_artifacts[d["dataset_id"]]
            for d in collection_meta["datasets"]
            if d["dataset_id"] in registered_artifacts
        ]
        if not member_artifacts:
            continue

        previous_collection = existing_pre_release.get(collection_meta["collection_id"])

        collection_kwargs: dict[str, Any] = {
            "key": f"pre-release/{collection_meta['name']}",
            "description": collection_meta["doi"],
            "reference": collection_meta["collection_id"],
            "reference_type": "CELLxGENE Collection ID",
        }
        if previous_collection is not None:
            collection_kwargs["revises"] = previous_collection
            logger.info(f"revising pre-release collection: {collection_meta['name']}")
        else:
            logger.info(f"creating pre-release collection: {collection_meta['name']}")

        collection_record = ln.Collection(member_artifacts, **collection_kwargs)
        collection_record.save()
        collection_record.ulabels.add(pre_release_label)

    ln.settings.creation.search_names = True


def _annotate_artifacts(
    cxg_datasets: list[dict[str, Any]],
    registered_ids: set[str],
    new_census_version: str,
    pre_release_label: Any = None,
) -> None:
    """Validate and curate registered artifacts."""
    logger.important("starting annotation of artifacts")

    # bulk query: fetch all organism records needed across all datasets at once
    all_organism_ontology_ids = {
        organism["ontology_term_id"]
        for ds in cxg_datasets
        if ds["dataset_id"] in registered_ids
        for organism in ds["organism"]
    }
    organism_lookup: dict[str, Any] = {
        o.ontology_id: o
        for o in bt.Organism.filter(ontology_id__in=all_organism_ontology_ids)
    }

    # pre-build schemas per unique organism — at most K requests (K = unique organisms)
    schema_cache: dict[str, Any] = {}
    for o in organism_lookup.values():
        name = "mouse" if o.name == "house mouse" else o.name
        if name in schema_cache:
            continue
        try:
            schema_cache[name] = ln.examples.cellxgene.create_cellxgene_schema(
                field_types="ontology_id", organism=name
            )
        except ObjectDoesNotExist:
            logger.warning(
                f"no schema for organism={name}: bt.Source not found. "
                f"run bt.Gene.add_source(organism='{name}') to fix."
            )
        except IndexError:
            logger.warning(
                f"no schema for organism={name}: IndexError while creating schema"
            )

    for idx, ds in enumerate(cxg_datasets):
        if ds["dataset_id"] not in registered_ids:
            continue
        if idx % 10 == 0:
            logger.info(f"annotating dataset {idx} of {len(cxg_datasets)}")

        if pre_release_label is not None:
            af = ln.Artifact.filter(
                ulabels=pre_release_label,
                key__endswith=f"{ds['dataset_id']}.h5ad",
            ).one_or_none()
        else:
            af = ln.Artifact.filter(
                Q(key__contains=ds["dataset_id"]) & Q(key__contains=new_census_version)
            ).one_or_none()

        if af is None:
            logger.warning(
                f"no artifact found for dataset_id={ds['dataset_id']}, skipping"
            )
            continue

        organism_ontology_ids = [
            organism["ontology_term_id"] for organism in ds["organism"]
        ]
        organism_records = [
            organism_lookup[oid]
            for oid in organism_ontology_ids
            if oid in organism_lookup
        ]

        if not organism_records:
            logger.warning(
                f"skipping dataset_id={ds['dataset_id']}: no organism records found"
            )
            continue
        first_organism = organism_records[0]
        organism_name = (
            "mouse" if first_organism.name == "house mouse" else first_organism.name
        )

        schema = schema_cache.get(organism_name)
        if schema is None:
            logger.warning(
                f"skipping dataset_id={ds['dataset_id']}: no schema for organism={organism_name}"
            )
            continue

        try:
            curator = ln.curators.AnnDataCurator(af, schema)
            curator.validate()
            curator.save_artifact(revises=af)
            logger.info(
                f"successfully validated and saved dataset_id={ds['dataset_id']}"
            )
        except ln.errors.ValidationError as e:
            error_msg = str(e)
            if "not validated in feature 'tissue_ontology_term_id'" in error_msg:
                logger.warning(
                    f"skipping dataset_id={ds['dataset_id']}: tissue_ontology_term_id not validated"
                )
            elif (
                "term not validated in feature 'self_reported_ethnicity_ontology_term_id' in slot 'obs'"
                in error_msg
            ):
                logger.warning(
                    f"skipping dataset_id={ds['dataset_id']}: self_reported_ethnicity_ontology_term_id not validated"
                )
            elif (
                "not validated in feature 'disease_ontology_term_id' in slot 'obs'"
                in error_msg
            ):
                logger.warning(
                    f"skipping dataset_id={ds['dataset_id']}: disease_ontology_term_id not validated"
                )
            elif "not in dataframe" in error_msg:
                logger.warning(
                    f"skipping dataset_id={ds['dataset_id']}: feature not in dataframe"
                )
            elif "no Organism found in source for the given" in error_msg:
                logger.warning(
                    f"skipping dataset_id={ds['dataset_id']}: ontology not in dataframe"
                )
            else:
                logger.error(
                    f"unhandled ValidationError for dataset_id={ds['dataset_id']}: {error_msg}"
                )


# ---------------------------------------------------------------------------
# Python API
# ---------------------------------------------------------------------------


@ln.flow("Rrq1bb328HH4")
def ingest_lts(
    new: str,
    previous: str | None = None,
    smoke: bool = False,
) -> None:
    """Register a full LTS census release into LaminDB."""
    logger.info(f"ingest-lts | new={new} | previous={previous} | smoke={smoke}")
    census_s3_path = f"s3://cellxgene-data-public/cell-census/{new}/h5ads"

    if smoke:
        ln.examples.cellxgene.save_cellxgene_defaults()
        logger.info("smoke mode: saved cellxgene defaults")

    cxg_datasets: list[dict[str, Any]] = get_datasets_from_cxg()  # type: ignore
    logger.info(f"found {len(cxg_datasets)} datasets from CellxGene")
    cxg_lookup: dict[str, dict[str, Any]] = {
        ds["dataset_id"]: ds for ds in cxg_datasets
    }  # type: ignore

    # 1. Register artifacts
    h5ad_paths = list(ln.UPath(census_s3_path).glob("*.h5ad"))
    logger.info(f"found {len(h5ad_paths)} h5ad paths in {census_s3_path}")
    if smoke:
        h5ad_paths = h5ad_paths[:2]

    registered_ids: set[str] = set()
    for h5ad_path in h5ad_paths:
        dataset_id = h5ad_path.stem
        registered_ids.add(dataset_id)
        artifact_previous = (
            ln.Artifact.filter(
                key__endswith=f"{dataset_id}.h5ad",
                key__contains=previous,
            ).one_or_none()
            if previous
            else None
        )
        kwargs: dict[str, Any] = {}
        if artifact_previous is not None:
            kwargs["revises"] = artifact_previous
            logger.info(f"revising existing artifact for dataset_id={dataset_id}")
        else:
            logger.info(f"registering new artifact for dataset_id={dataset_id}")

        artifact = ln.Artifact(h5ad_path, **kwargs)
        artifact.version_tag = new
        if dataset_id in cxg_lookup:
            artifact.description = cxg_lookup[dataset_id]["title"]
            artifact.n_observations = cxg_lookup[dataset_id]["cell_count"]
        artifact.save()

    new_afs = ln.Artifact.filter(key__contains=new)
    logger.info(f"registered {len(registered_ids)} artifacts for census version {new}")

    if not smoke:
        # 2. Register collections
        logger.info("registering top-level cellxgene-census collection")
        collection = ln.Collection(
            new_afs,
            key="cellxgene-census",
            revises=ln.Collection.filter(
                key="cellxgene-census", version_tag=previous
            ).one(),
        )
        collection.version_tag = new
        collection.save()
        logger.info("saved top-level collection")

        cxg_collections: list[dict[str, Any]] = get_collections_from_cxg()  # type: ignore[assignment]
        logger.info(f"found {len(cxg_collections)} CellxGene collections")

        ln.settings.creation.search_names = False
        for collection_meta in cxg_collections:
            keys = [
                f"cell-census/{new}/h5ads/{dataset['dataset_id']}.h5ad"
                for dataset in collection_meta["datasets"]
            ]
            collection_artifacts = new_afs.filter(key__in=keys)
            if collection_artifacts.count() > 0:
                previous_collection = ln.Collection.filter(
                    reference=collection_meta["collection_id"],
                    version_tag=previous,
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
                collection_record.version_tag = new
                collection_record.save()
            else:
                logger.warning(
                    f"no matching artifacts for collection: {collection_meta['name']} (id={collection_meta['collection_id']}), skipping"
                )
        ln.settings.creation.search_names = True

        # 3. Register soma store
        logger.info("registering soma store")
        soma_path = f"s3://cellxgene-data-public/cell-census/{new}/soma"
        previous_soma = ln.Artifact.filter(description=f"Census {previous}").one()
        new_soma_af = ln.Artifact(
            soma_path,
            description=f"Census {new}",
            revises=previous_soma,
        )
        new_soma_af.version_tag = new
        new_soma_af.save()
        logger.info(f"saved soma artifact: {new_soma_af}")

    # 4. Annotate (skipped in smoke mode)
    if not smoke:
        _annotate_artifacts(cxg_datasets, registered_ids, new_census_version=new)


@ln.flow("Rrq1bb328HH4")
def ingest_pre_release(
    new: str,
    smoke: bool = False,
) -> None:
    """Register datasets from the latest weekly census not yet in LaminDB."""
    import cellxgene_census as cxc

    logger.info(f"ingest-pre-release | new={new} | smoke={smoke}")

    if smoke:
        ln.examples.cellxgene.save_cellxgene_defaults()
        logger.info("smoke mode: saved cellxgene defaults")

    cxg_datasets: list[dict[str, Any]] = get_datasets_from_cxg()  # type: ignore
    logger.info(f"found {len(cxg_datasets)} datasets from CellxGene")
    cxg_lookup: dict[str, dict[str, Any]] = {
        ds["dataset_id"]: ds for ds in cxg_datasets
    }  # type: ignore

    # get or create pre-release label
    pre_release_label = ln.ULabel.filter(name="pre-release").one_or_none()
    if pre_release_label is None:
        pre_release_label = ln.ULabel(name="pre-release").save()
    logger.info("pre-release ULabel ready")

    # 1. Cleanup broken artifacts from previous runs
    cleanup_broken_pre_release_artifacts(pre_release_label)

    # 2. Register new artifacts
    registered_ids, registered_artifacts = register_pre_release_artifacts(
        cxg_lookup=cxg_lookup,
        pre_release_label=pre_release_label,
        cxc=cxc,
        smoke=smoke,
    )
    logger.info(f"registered {len(registered_ids)} pre-release artifacts")

    # 3. Register collections not in LTS
    if registered_artifacts and not smoke:
        register_pre_release_collections(
            cxg_collections=get_collections_from_cxg(),  # type: ignore[arg-type]
            registered_artifacts=registered_artifacts,
            pre_release_label=pre_release_label,
        )

    # 4. Annotate (skipped in smoke mode)
    if not smoke:
        _annotate_artifacts(
            cxg_datasets,
            registered_ids,
            new_census_version=new,
            pre_release_label=pre_release_label,
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.group()
def cxg_lamin() -> None:
    """CellxGene LaminDB CLI."""


@cxg_lamin.command("ingest-lts")
@click.option("--new", required=True, help="New census version")
@click.option("--previous", default=None, help="Previous census version")
@click.option("--smoke", is_flag=True, help="Limit to 2 datasets, skip collections")
def ingest_lts_cmd(new: str, previous: str | None, smoke: bool) -> None:
    """Register a full LTS census release."""
    ingest_lts(new=new, previous=previous, smoke=smoke)


@cxg_lamin.command("ingest-pre-release")
@click.option("--new", required=True, help="New census version")
@click.option("--smoke", is_flag=True, help="Limit to 2 datasets, skip collections")
def ingest_pre_release_cmd(new: str, smoke: bool) -> None:
    """Register pre-release datasets from the latest weekly census."""
    ingest_pre_release(new=new, smoke=smoke)


def main() -> None:
    """Entry point for the cxg-lamin CLI."""
    cxg_lamin()
