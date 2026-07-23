"""Tests for register-annotate-new-release.py"""

import importlib.util
import pathlib
import sys
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# Keep a reference to the lamindb mock so tests can patch attributes on it.
# patch.dict restores sys.modules after the block, so we must hold the mock
# ourselves to reach the same object the module's `ln` points to.
# ---------------------------------------------------------------------------
_lamindb_mock = MagicMock()
_mocks = {
    "lamindb": _lamindb_mock,
    "bionty": MagicMock(),
    "lamin_utils": MagicMock(),
    "django": MagicMock(),
    "django.db": MagicMock(),
    "django.db.models": MagicMock(),
    "django.core": MagicMock(),
    "django.core.exceptions": MagicMock(),
    "cellxgene_lamin": MagicMock(),
    "cellxgene_lamin.dev": MagicMock(),
    "cellxgene_lamin.dev._cxg_rest_api": MagicMock(),
}

with (
    patch.dict("sys.modules", _mocks),
    patch("sys.argv", ["script", "--new", "2025-07-01"]),
):
    _spec = importlib.util.spec_from_file_location(
        "register_annotate",
        pathlib.Path(__file__).parent.parent
        / "cellxgene_lamin/_register_annotate_new_release.py",
    )
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    register_pre_release_artifacts = _mod.register_pre_release_artifacts
    cleanup_broken_pre_release_artifacts = _mod.cleanup_broken_pre_release_artifacts
    register_pre_release_collections = _mod.register_pre_release_collections
    ingest_lts = _mod.ingest_lts


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def make_ds(dataset_id="ds-001", title="Test Dataset", cell_count=100):
    return {
        "dataset_id": dataset_id,
        "title": title,
        "cell_count": cell_count,
        "organism": [{"ontology_term_id": "NCBITaxon:10090"}],
    }


# ---------------------------------------------------------------------------
# Test 1: dataset already in LaminDB is skipped
# ---------------------------------------------------------------------------


def test_skips_already_registered_dataset():
    """If an artifact with the dataset's key already exists, it must be skipped."""
    cxg_lookup = {"ds-001": make_ds("ds-001")}
    pre_release_label = MagicMock()
    mock_cxc = MagicMock()

    with patch.object(_lamindb_mock, "Artifact") as mock_artifact:
        # Simulate: artifact already exists in LaminDB
        mock_artifact.filter.return_value.one_or_none.return_value = MagicMock()

        registered_ids, registered_artifacts = register_pre_release_artifacts(
            cxg_lookup=cxg_lookup,
            pre_release_label=pre_release_label,
            cxc=mock_cxc,
        )

    # dataset must not appear in the registered sets
    assert "ds-001" not in registered_ids
    assert "ds-001" not in registered_artifacts

    # URI resolution should never have been attempted
    mock_cxc.get_source_h5ad_uri.assert_not_called()


# ---------------------------------------------------------------------------
# Test 2: dataset not in LaminDB, URI resolves, file accessible → registered
# ---------------------------------------------------------------------------


def test_registers_new_dataset():
    """If a dataset is not in LaminDB and the file is accessible, it gets registered."""
    cxg_lookup = {"ds-001": make_ds("ds-001")}
    pre_release_label = MagicMock()
    mock_cxc = MagicMock()
    mock_cxc.get_source_h5ad_uri.return_value = {"uri": "s3://bucket/ds-001.h5ad"}

    with patch.object(_lamindb_mock, "Artifact") as mock_artifact_class:
        # dataset not in LaminDB
        mock_artifact_class.filter.return_value.one_or_none.return_value = None

        # the artifact instance the constructor returns
        mock_artifact_instance = MagicMock()
        mock_artifact_class.return_value = mock_artifact_instance

        registered_ids, registered_artifacts = register_pre_release_artifacts(
            cxg_lookup=cxg_lookup,
            pre_release_label=pre_release_label,
            cxc=mock_cxc,
        )

    assert "ds-001" in registered_ids
    assert "ds-001" in registered_artifacts
    mock_artifact_instance.save.assert_called_once()
    mock_artifact_instance.ulabels.add.assert_called_once_with(pre_release_label)


# ---------------------------------------------------------------------------
# Test 3: URI resolution fails → dataset skipped
# ---------------------------------------------------------------------------


def test_skips_dataset_when_uri_resolution_fails():
    """If get_source_h5ad_uri raises, the dataset must be skipped."""
    cxg_lookup = {"ds-001": make_ds("ds-001")}
    pre_release_label = MagicMock()
    mock_cxc = MagicMock()
    mock_cxc.get_source_h5ad_uri.side_effect = Exception("dataset not in census")

    with patch.object(_lamindb_mock, "Artifact") as mock_artifact_class:
        # dataset not in LaminDB
        mock_artifact_class.filter.return_value.one_or_none.return_value = None
        mock_artifact_instance = MagicMock()
        mock_artifact_class.return_value = mock_artifact_instance

        registered_ids, registered_artifacts = register_pre_release_artifacts(
            cxg_lookup=cxg_lookup,
            pre_release_label=pre_release_label,
            cxc=mock_cxc,
        )

    assert "ds-001" not in registered_ids
    mock_artifact_instance.save.assert_not_called()


# ---------------------------------------------------------------------------
# Test 4: artifact.open() raises → dataset skipped before save
# ---------------------------------------------------------------------------


def test_skips_broken_dataset_before_registering():
    """If artifact.open() raises, the dataset must be skipped and never saved."""
    cxg_lookup = {"ds-001": make_ds("ds-001")}
    pre_release_label = MagicMock()
    mock_cxc = MagicMock()
    mock_cxc.get_source_h5ad_uri.return_value = {"uri": "s3://bucket/ds-001.h5ad"}

    with patch.object(_lamindb_mock, "Artifact") as mock_artifact_class:
        # dataset not in LaminDB
        mock_artifact_class.filter.return_value.one_or_none.return_value = None

        # artifact instance whose open() raises (broken file)
        mock_artifact_instance = MagicMock()
        mock_artifact_instance.open.side_effect = Exception("file not accessible")
        mock_artifact_class.return_value = mock_artifact_instance

        registered_ids, registered_artifacts = register_pre_release_artifacts(
            cxg_lookup=cxg_lookup,
            pre_release_label=pre_release_label,
            cxc=mock_cxc,
        )

    assert "ds-001" not in registered_ids
    mock_artifact_instance.save.assert_not_called()
    mock_artifact_instance.ulabels.add.assert_not_called()


# ---------------------------------------------------------------------------
# Test 5: existing pre-release artifact is broken → cleanup deletes it
# ---------------------------------------------------------------------------


def test_cleanup_deletes_broken_pre_release_artifact():
    """Broken pre-release artifacts from previous runs must be permanently deleted."""
    pre_release_label = MagicMock()

    broken_artifact = MagicMock()
    broken_artifact.open.side_effect = Exception("file not accessible")

    ok_artifact = MagicMock()  # open() does not raise → should NOT be deleted

    with patch.object(_lamindb_mock, "Artifact") as mock_artifact_class:
        mock_artifact_class.filter.return_value.count.return_value = 2
        mock_artifact_class.filter.return_value.__iter__ = lambda s: iter(
            [broken_artifact, ok_artifact]
        )

        cleanup_broken_pre_release_artifacts(pre_release_label)

    broken_artifact.delete.assert_called_once_with(permanent=True)
    ok_artifact.delete.assert_not_called()


# ---------------------------------------------------------------------------
# Test 6: collection already in LTS → skipped, no Collection saved
# ---------------------------------------------------------------------------


def make_collection_meta(collection_id="col-001", name="My Collection", datasets=None):
    return {
        "collection_id": collection_id,
        "name": name,
        "doi": "10.1234/test",
        "datasets": datasets or [{"dataset_id": "ds-001"}],
    }


def test_skips_collection_already_in_lts():
    """Collections that already have a version-tagged record in LaminDB must be skipped."""
    pre_release_label = MagicMock()
    registered_artifacts = {"ds-001": MagicMock()}
    cxg_collections = [make_collection_meta("col-001")]

    with (
        patch.object(_lamindb_mock, "Collection") as mock_collection_class,
        patch.object(_lamindb_mock, "settings"),
    ):
        # LTS check returns True → collection already in LTS
        mock_collection_class.filter.return_value.exists.return_value = True

        register_pre_release_collections(
            cxg_collections=cxg_collections,
            registered_artifacts=registered_artifacts,
            pre_release_label=pre_release_label,
        )

    # no Collection should have been created or saved
    mock_collection_class.assert_not_called()
    mock_collection_class.return_value.save.assert_not_called()


# ---------------------------------------------------------------------------
# Test 7: brand-new collection → registered with pre-release/ key and label
# ---------------------------------------------------------------------------


def test_registers_new_pre_release_collection():
    """A collection not in LTS and not previously registered must be created and labeled."""
    pre_release_label = MagicMock()
    ds_artifact = MagicMock()
    registered_artifacts = {"ds-001": ds_artifact}
    cxg_collections = [
        make_collection_meta("col-001", "My Collection", [{"dataset_id": "ds-001"}])
    ]

    with (
        patch.object(_lamindb_mock, "Collection") as mock_collection_class,
        patch.object(_lamindb_mock, "settings"),
    ):
        # not in LTS
        mock_collection_class.filter.return_value.exists.return_value = False
        # no previous pre-release collection
        mock_collection_class.filter.return_value.one_or_none.return_value = None

        mock_collection_instance = MagicMock()
        mock_collection_class.return_value = mock_collection_instance

        register_pre_release_collections(
            cxg_collections=cxg_collections,
            registered_artifacts=registered_artifacts,
            pre_release_label=pre_release_label,
        )

    # Collection created with the right key
    mock_collection_class.assert_called_once_with(
        [ds_artifact],
        key="pre-release/My Collection",
        description="10.1234/test",
        reference="col-001",
        reference_type="CELLxGENE Collection ID",
    )
    mock_collection_instance.save.assert_called_once()
    mock_collection_instance.ulabels.add.assert_called_once_with(pre_release_label)


# ===========================================================================
# LTS path  (ingest_lts)
# ===========================================================================


def _make_h5ad_path(dataset_id: str):
    """Return a MagicMock that behaves like a UPath to '<dataset_id>.h5ad'."""
    p = MagicMock()
    p.stem = dataset_id
    return p


# ---------------------------------------------------------------------------
# Test 8: new LTS dataset (no previous version) → registered fresh
# ---------------------------------------------------------------------------


def test_lts_registers_new_artifact():
    """LTS: a dataset absent from the previous release is saved without revises."""
    mock_path = _make_h5ad_path("ds-001")

    with (
        patch.object(_mod, "get_datasets_from_cxg", return_value=[make_ds("ds-001")]),
        patch.object(_mod, "get_collections_from_cxg", return_value=[]),
        patch.object(_mod, "_annotate_artifacts"),  # skip annotation
        patch.object(_lamindb_mock, "Artifact") as mock_artifact_class,
        patch.object(_lamindb_mock, "UPath") as mock_upath,
        patch.object(_lamindb_mock, "Collection"),
    ):
        mock_upath.return_value.glob.return_value = [mock_path]
        # no previous artifact
        mock_artifact_class.filter.return_value.one_or_none.return_value = None

        mock_instance = MagicMock()
        mock_artifact_class.return_value = mock_instance

        ingest_lts(new="2025-07-01", previous="2025-01-30", smoke=True)

    # Artifact must have been constructed with the h5ad path (no revises)
    mock_artifact_class.assert_called_once_with(mock_path)
    assert mock_instance.version_tag == "2025-07-01"
    mock_instance.save.assert_called_once()


# ---------------------------------------------------------------------------
# Test 9: LTS dataset also in previous release → registered with revises=…
# ---------------------------------------------------------------------------


def test_lts_revises_previous_artifact():
    """LTS: a dataset present in the previous release is registered with revises=prev."""
    mock_path = _make_h5ad_path("ds-001")
    prev_artifact = MagicMock(name="prev_artifact")

    with (
        patch.object(_mod, "get_datasets_from_cxg", return_value=[make_ds("ds-001")]),
        patch.object(_mod, "get_collections_from_cxg", return_value=[]),
        patch.object(_mod, "_annotate_artifacts"),
        patch.object(_lamindb_mock, "Artifact") as mock_artifact_class,
        patch.object(_lamindb_mock, "UPath") as mock_upath,
        patch.object(_lamindb_mock, "Collection"),
    ):
        mock_upath.return_value.glob.return_value = [mock_path]
        # previous artifact found
        mock_artifact_class.filter.return_value.one_or_none.return_value = prev_artifact

        mock_instance = MagicMock()
        mock_artifact_class.return_value = mock_instance

        ingest_lts(new="2025-07-01", previous="2025-01-30", smoke=True)

    # Artifact must have been constructed with revises=prev_artifact
    mock_artifact_class.assert_called_once_with(mock_path, revises=prev_artifact)
    assert mock_instance.version_tag == "2025-07-01"
    mock_instance.save.assert_called_once()


# ---------------------------------------------------------------------------
# Test 10: smoke mode caps processing to the first 2 h5ad paths
# ---------------------------------------------------------------------------


def test_lts_smoke_limits_to_two_paths():
    """LTS smoke mode: only the first 2 paths from S3 are processed."""
    paths = [_make_h5ad_path(f"ds-{i:03d}") for i in range(5)]

    with (
        patch.object(_mod, "get_datasets_from_cxg", return_value=[]),
        patch.object(_mod, "_annotate_artifacts"),
        patch.object(_lamindb_mock, "Artifact") as mock_artifact_class,
        patch.object(_lamindb_mock, "UPath") as mock_upath,
    ):
        mock_upath.return_value.glob.return_value = paths
        mock_artifact_class.filter.return_value.one_or_none.return_value = None

        ingest_lts(new="2025-07-01", smoke=True)

    # With smoke=True the list is sliced to [:2] before the loop
    assert mock_artifact_class.call_count == 2


# ---------------------------------------------------------------------------
# Test 11: LTS no-smoke → top-level collection is created and version-tagged
# ---------------------------------------------------------------------------


def test_lts_non_smoke_creates_top_level_collection():
    """LTS full run: a top-level 'cellxgene-census' collection is saved with version_tag."""
    mock_path = _make_h5ad_path("ds-001")

    with (
        patch.object(_mod, "get_datasets_from_cxg", return_value=[make_ds("ds-001")]),
        patch.object(_mod, "get_collections_from_cxg", return_value=[]),
        patch.object(_mod, "_annotate_artifacts"),
        patch.object(_lamindb_mock, "Artifact") as mock_artifact_class,
        patch.object(_lamindb_mock, "UPath") as mock_upath,
        patch.object(_lamindb_mock, "Collection") as mock_collection_class,
    ):
        mock_upath.return_value.glob.return_value = [mock_path]
        mock_artifact_class.filter.return_value.one_or_none.return_value = None

        mock_collection_instance = MagicMock()
        mock_collection_class.return_value = mock_collection_instance

        ingest_lts(new="2025-07-01", previous="2025-01-30", smoke=False)

    # A top-level Collection(new_afs, key="cellxgene-census", ...) must be created
    assert mock_collection_class.call_count >= 1
    first_call_kwargs = mock_collection_class.call_args_list[0].kwargs
    assert first_call_kwargs.get("key") == "cellxgene-census"
    assert mock_collection_instance.version_tag == "2025-07-01"
    mock_collection_instance.save.assert_called()
