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
    """Delete pre-release artifacts only when confirmed missing from storage.

    Uses FileNotFoundError and botocore ClientError with NoSuchKey/NoSuchBucket
    codes to confirm the file is gone. Transient errors (timeouts, permission
    errors) are logged as warnings and do not trigger deletion.
    """
    try:
        from botocore.exceptions import ClientError as BotocoreClientError

        _missing_codes = {"NoSuchKey", "NoSuchBucket", "404"}

        def _is_missing(e: Exception) -> bool:
            if isinstance(e, FileNotFoundError):
                return True
            if isinstance(e, BotocoreClientError):
                return e.response["Error"]["Code"] in _missing_codes
            return False

    except ImportError:
        # botocore not available — fall back to FileNotFoundError only
        def _is_missing(e: Exception) -> bool:  # type: ignore[misc]
            return isinstance(e, FileNotFoundError)

    existing = ln.Artifact.filter(ulabels=pre_release_label)
    logger.info(
        f"checking {existing.count()} existing pre-release artifacts for breakage"
    )
    for af in existing:
        try:
            af.open()
        except Exception as e:
            if _is_missing(e):
                logger.warning(f"confirmed missing, deleting | key={af.key} | {e}")
                af.delete(permanent=True)
            else:
                logger.warning(f"transient error, skipping | key={af.key} | {e}")


def _build_census_uri_map(cxc: Any) -> dict[str, str]:
    """Return a mapping of census dataset_id → full S3 h5ad URI."""
    with cxc.open_soma(census_version="latest") as census:
        h5ad_base = census.uri.replace("/soma", "/h5ads").rstrip("/")
        census_df = (
            census["census_info"]["datasets"]
            .read()
            .concat()
            .to_pandas()[["dataset_id", "dataset_h5ad_path"]]
        )
    return {
        row["dataset_id"]: f"{h5ad_base}/{row['dataset_h5ad_path']}"
        for _, row in census_df.iterrows()
    }


def register_pre_release_artifacts(
    cxg_lookup: dict[str, dict[str, Any]],
    pre_release_label: Any,
    cxc: Any,
    n_datasets: int = 5,
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
    census_uri_map = _build_census_uri_map(cxc)
    logger.info(f"found {len(census_uri_map)} datasets in census")

    cxg_lookup_filtered = {
        dataset_id: ds
        for dataset_id, ds in cxg_lookup.items()
        if dataset_id not in existing_dataset_ids and dataset_id in census_uri_map
    }
    logger.info(f"{len(cxg_lookup_filtered)} datasets to register")

    items = list(cxg_lookup_filtered.items())
    if n_datasets >= 0:
        items = items[:n_datasets]

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


def register_portal_only_artifacts(
    cxg_lookup: dict[str, dict[str, Any]],
    pre_release_label: Any,
    portal_only_label: Any,
    census_uri_map: dict[str, str],
    n_datasets: int = 5,
) -> tuple[set[str], dict[str, Any]]:
    """Register datasets that are on the CellxGENE portal but not in Census.

    These datasets have no stable S3 URI in the Census, so they are registered
    via their HTTPS download URL from the portal (``assets[].url``). They are
    tagged with both the ``pre-release`` and ``portal-only`` ULabels.

    Deduplication uses the ``dataset_version_id`` embedded in the download URL
    (the filename stem), because that is what appears in the artifact key —
    not the ``dataset_id``.
    """
    registered_ids: set[str] = set()
    registered_artifacts: dict[str, Any] = {}

    existing_version_ids = {
        key.rsplit("/", 1)[-1].removesuffix(".h5ad")
        for key in ln.Artifact.filter(ulabels=portal_only_label).values_list(
            "key", flat=True
        )
    }
    logger.info(
        f"found {len(existing_version_ids)} already-registered portal-only datasets"
    )

    candidates: list[tuple[str, str, str, dict[str, Any]]] = []
    for dataset_id, ds in cxg_lookup.items():
        if dataset_id in census_uri_map:
            continue
        h5ad_url = next(
            (a["url"] for a in ds.get("assets", []) if a.get("filetype") == "H5AD"),
            None,
        )
        if h5ad_url is None:
            continue
        dataset_version_id = h5ad_url.rsplit("/", 1)[-1].removesuffix(".h5ad")
        if dataset_version_id in existing_version_ids:
            continue
        candidates.append((dataset_id, dataset_version_id, h5ad_url, ds))

    logger.info(
        f"found {len(candidates)} portal-only datasets to register (not in Census, not yet in LaminDB)"
    )

    if n_datasets >= 0:
        candidates = candidates[:n_datasets]

    for dataset_id, _version_id, h5ad_url, ds in candidates:
        artifact = ln.Artifact(h5ad_url, description=ds["title"])
        artifact.n_observations = ds["cell_count"]
        artifact.save()
        artifact.ulabels.add(pre_release_label)
        artifact.ulabels.add(portal_only_label)

        registered_ids.add(dataset_id)
        registered_artifacts[dataset_id] = artifact
        logger.info(f"registered portal-only | dataset_id={dataset_id}")

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

    # build a lookup of all existing pre-release artifacts by dataset_id
    all_pre_release_artifacts: dict[str, Any] = {
        af.key.rsplit("/", 1)[-1].removesuffix(".h5ad"): af
        for af in ln.Artifact.filter(ulabels=pre_release_label)
    }
    # merge with newly registered ones (newly registered take precedence)
    all_pre_release_artifacts.update(registered_artifacts)

    ln.settings.creation.search_names = False
    for collection_meta in cxg_collections:
        if collection_meta["collection_id"] in lts_collection_ids:
            logger.info(
                f"collection already in LTS, skipping | id={collection_meta['collection_id']}"
            )
            continue

        member_artifacts = [
            all_pre_release_artifacts[d["dataset_id"]]
            for d in collection_meta["datasets"]
            if d["dataset_id"] in all_pre_release_artifacts
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
    artifact_lookup: dict[str, Any] | None = None,
) -> None:
    """Validate and curate registered artifacts.

    Args:
        cxg_datasets: Dataset metadata records from the CellxGENE portal API.
        registered_ids: Dataset IDs to annotate.
        new_census_version: Census version date used to look up LTS artifacts.
        pre_release_label: ULabel used to look up pre-release artifacts. When
            omitted, artifacts are matched by census version instead.
        artifact_lookup: Optional direct mapping of ``dataset_id`` to artifact.
            When provided, it is checked first before falling back to the
            key-based query. Required for portal-only datasets whose key stem
            is the ``dataset_version_id``, not the ``dataset_id``.
    """
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

        if artifact_lookup is not None and ds["dataset_id"] in artifact_lookup:
            af = artifact_lookup[ds["dataset_id"]]
        elif pre_release_label is not None:
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
        except Exception as e:
            # transient errors (network timeouts, socket errors) — log and skip,
            # the dataset will be retried on the next weekly run via the
            # unannotated_ids recovery path
            logger.warning(
                f"transient error annotating dataset_id={ds['dataset_id']}, will retry next run | {type(e).__name__}: {e}"
            )


# ---------------------------------------------------------------------------
# Python API
# ---------------------------------------------------------------------------


@ln.flow("Rrq1bb328HH4")
def ingest_lts(
    new: str,
    previous: str = "",
    n_datasets: int = 5,
) -> None:
    """Register a full LTS census release into LaminDB.

    Args:
        new: New census version date (e.g. '2025-11-08').
        previous: Previous census version date, used to revise existing artifacts.
        n_datasets: Max number of datasets to ingest. Pass -1 to ingest all
            (including collections and soma store). Defaults to 5.
    """
    is_full_run = n_datasets < 0
    logger.info(
        f"ingest-lts | new={new} | previous={previous} | n_datasets={n_datasets}"
    )
    census_s3_path = f"s3://cellxgene-data-public/cell-census/{new}/h5ads"

    if not is_full_run:
        ln.examples.cellxgene.save_cellxgene_defaults()

    cxg_datasets: list[dict[str, Any]] = get_datasets_from_cxg()  # type: ignore
    logger.info(f"found {len(cxg_datasets)} datasets from CellxGene")
    cxg_lookup: dict[str, dict[str, Any]] = {
        ds["dataset_id"]: ds for ds in cxg_datasets
    }  # type: ignore

    # 1. Register artifacts
    h5ad_paths = list(ln.UPath(census_s3_path).glob("*.h5ad"))
    logger.info(f"found {len(h5ad_paths)} h5ad paths in {census_s3_path}")
    if not is_full_run:
        h5ad_paths = h5ad_paths[:n_datasets]

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

    # strip pre-release label from any artifacts now promoted to LTS
    pre_release_label = ln.ULabel.filter(name="pre-release").one_or_none()
    if pre_release_label is not None:
        lts_dataset_ids = {
            key.rsplit("/", 1)[-1].removesuffix(".h5ad")
            for key in ln.Artifact.filter(version_tag=new).values_list("key", flat=True)
        }
        count = 0
        for dataset_id in lts_dataset_ids:
            af = ln.Artifact.filter(
                ulabels=pre_release_label,
                key__endswith=f"{dataset_id}.h5ad",
            ).one_or_none()
            if af is not None:
                af.ulabels.remove(pre_release_label)
                count += 1
        if count:
            logger.info(
                f"stripped pre-release label from {count} artifacts now in LTS {new}"
            )

    if is_full_run:
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

    # 4. Annotate
    _annotate_artifacts(cxg_datasets, registered_ids, new_census_version=new)


@ln.flow("Rrq1bb328HH4")
def ingest_pre_release(
    new: str,
    n_datasets: int = 5,
) -> None:
    """Register datasets from the latest weekly census not yet in LaminDB.

    Args:
        new: LTS census version date used to identify already-ingested datasets.
        n_datasets: Max number of new pre-release datasets to register. Pass -1
            to register all (including collections). Defaults to 5.
    """
    import cellxgene_census as cxc

    is_full_run = n_datasets < 0
    logger.info(f"ingest-pre-release | new={new} | n_datasets={n_datasets}")

    if not is_full_run:
        ln.examples.cellxgene.save_cellxgene_defaults()

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
        n_datasets=n_datasets,
    )
    logger.info(f"registered {len(registered_ids)} pre-release artifacts")

    # 3. Register collections not in LTS (full run only)
    if registered_artifacts and is_full_run:
        register_pre_release_collections(
            cxg_collections=get_collections_from_cxg(),  # type: ignore[arg-type]
            registered_artifacts=registered_artifacts,
            pre_release_label=pre_release_label,
        )

    # 4. Annotate — include any pre-release artifacts that were previously
    # registered but not annotated (e.g. a prior run stopped mid-way)
    unannotated_ids = {
        af.key.rsplit("/", 1)[-1].removesuffix(".h5ad")
        for af in ln.Artifact.filter(ulabels=pre_release_label)
        if af.cell_types.count() == 0 and af.tissues.count() == 0
    }
    ids_to_annotate = registered_ids | unannotated_ids
    if unannotated_ids:
        logger.info(
            f"found {len(unannotated_ids)} pre-release artifacts missing annotation — will annotate now"
        )
    _annotate_artifacts(
        cxg_datasets,
        ids_to_annotate,
        new_census_version=new,
        pre_release_label=pre_release_label,
    )


@ln.flow("Rrq1bb328HH4")
def ingest_portal_only(
    new: str,
    n_datasets: int = 5,
) -> None:
    """Register datasets that are on the CellxGENE portal but not in Census.

    Datasets published on the portal but not yet included in the Census
    (for example spatial datasets or very new uploads) have no stable S3 URI.
    This function registers them via their HTTPS download URL and tags them
    with both ``pre-release`` and ``portal-only`` ULabels.

    Args:
        new: Current LTS census version date (e.g. ``'2025-11-08'``). Used
            only to identify the Census exclusion list — datasets already in
            the Census are skipped.
        n_datasets: Max number of new portal-only datasets to register. Pass
            ``-1`` to register all. Defaults to 5 for a quick test run.
    """
    import cellxgene_census as cxc

    is_full_run = n_datasets < 0
    logger.info(f"ingest-portal-only | new={new} | n_datasets={n_datasets}")

    if not is_full_run:
        ln.examples.cellxgene.save_cellxgene_defaults()

    cxg_datasets: list[dict[str, Any]] = get_datasets_from_cxg()  # type: ignore
    logger.info(f"found {len(cxg_datasets)} datasets from CellxGene portal")
    cxg_lookup: dict[str, dict[str, Any]] = {
        ds["dataset_id"]: ds for ds in cxg_datasets
    }

    pre_release_label = ln.ULabel.filter(name="pre-release").one_or_none()
    if pre_release_label is None:
        pre_release_label = ln.ULabel(name="pre-release").save()

    portal_only_label = ln.ULabel.filter(name="portal-only").one_or_none()
    if portal_only_label is None:
        portal_only_label = ln.ULabel(name="portal-only").save()
    logger.info("ULabels ready")

    census_uri_map = _build_census_uri_map(cxc)
    logger.info(f"Census has {len(census_uri_map)} datasets (will be excluded)")

    registered_ids, registered_artifacts = register_portal_only_artifacts(
        cxg_lookup=cxg_lookup,
        pre_release_label=pre_release_label,
        portal_only_label=portal_only_label,
        census_uri_map=census_uri_map,
        n_datasets=n_datasets,
    )
    logger.info(f"registered {len(registered_ids)} new portal-only artifacts")

    # The key stem is the dataset_version_id, so rebuild dataset_id → artifact
    # for any portal-only artifacts that were registered but not annotated.
    version_id_to_dataset_id: dict[str, str] = {}
    for dataset_id, ds in cxg_lookup.items():
        for asset in ds.get("assets", []):
            if asset.get("filetype") == "H5AD":
                vid = asset["url"].rsplit("/", 1)[-1].removesuffix(".h5ad")
                version_id_to_dataset_id[vid] = dataset_id
                break

    unannotated_artifact_lookup: dict[str, Any] = {}
    for af in ln.Artifact.filter(ulabels=portal_only_label):
        if af.cell_types.count() == 0 and af.tissues.count() == 0:
            vid = af.key.rsplit("/", 1)[-1].removesuffix(".h5ad")
            did = version_id_to_dataset_id.get(vid)
            if did:
                unannotated_artifact_lookup[did] = af

    if unannotated_artifact_lookup:
        logger.info(
            f"found {len(unannotated_artifact_lookup)} portal-only artifacts missing annotation — will annotate now"
        )

    combined_lookup = {**unannotated_artifact_lookup, **registered_artifacts}
    ids_to_annotate = set(combined_lookup.keys())

    _annotate_artifacts(
        cxg_datasets,
        ids_to_annotate,
        new_census_version=new,
        pre_release_label=portal_only_label,
        artifact_lookup=combined_lookup,
    )


# ---------------------------------------------------------------------------
# LaminHub automation entry point
# ---------------------------------------------------------------------------


@ln.flow()
def lamin_executable_ingest_pre_release(
    new: str = "2025-11-08",
    n_datasets: int = -1,
) -> None:
    """Weekly pre-release ingestion job — safe to re-run (idempotent).

    Registers datasets from the latest CellxGENE Census build that are not yet
    part of the LTS release identified by `new`. Annotates each new artifact
    with ontology-backed labels. Cleans up broken pre-release artifacts before
    registering new ones.

    Args:
        new: Current LTS census version date (e.g. '2025-11-08'). Update this
            when a new LTS is released.
        n_datasets: Max datasets to register per run. -1 registers all (default
            for weekly job). Pass a small number (e.g. 5) for a test run.

    Idempotency guarantees:
        - Already-registered datasets are detected via key lookup and skipped.
        - Broken pre-release artifacts are deleted before registering new ones.
        - Running twice produces no duplicates and deletes no healthy artifacts.
    """
    ingest_pre_release(new=new, n_datasets=n_datasets)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.group()
def cxg_lamin() -> None:
    """CellxGene LaminDB CLI."""


@cxg_lamin.command("ingest-lts")
@click.option("--new", required=True, help="New census version")
@click.option("--previous", default="", help="Previous census version")
@click.option(
    "--n-datasets",
    default=5,
    show_default=True,
    help="Max datasets to ingest with annotation; -1 ingests all (including collections and soma store).",
)
def ingest_lts_cmd(new: str, previous: str, n_datasets: int) -> None:
    """Register a full LTS census release."""
    ingest_lts(new=new, previous=previous, n_datasets=n_datasets)


@cxg_lamin.command("ingest-pre-release")
@click.option("--new", required=True, help="New census version")
@click.option(
    "--n-datasets",
    default=5,
    show_default=True,
    help="Max datasets to ingest with annotation; -1 ingests all (including collections).",
)
def ingest_pre_release_cmd(new: str, n_datasets: int) -> None:
    """Register pre-release datasets from the latest weekly census."""
    ingest_pre_release(new=new, n_datasets=n_datasets)


@cxg_lamin.command("ingest-portal-only")
@click.option("--new", required=True, help="Current LTS census version date")
@click.option(
    "--n-datasets",
    default=5,
    show_default=True,
    help="Max portal-only datasets to ingest; -1 registers all.",
)
def ingest_portal_only_cmd(new: str, n_datasets: int) -> None:
    """Register portal datasets not yet in Census."""
    ingest_portal_only(new=new, n_datasets=n_datasets)


def main() -> None:
    """Entry point for the cxg-lamin CLI."""
    cxg_lamin()
