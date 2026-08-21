"""Integration smoke tests — run against a fresh local lamindb instance.

A fresh local SQLite instance is used so the tests are fully isolated and
don't pollute any shared instance. The module is imported inside each test
function (not at module level) so Django is fully initialized by
ln.setup.init before lamindb models are loaded.
"""

import shutil
import sys
from pathlib import Path

import lamindb as ln
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

LTS_NEW = "2025-11-08"
LTS_PREVIOUS = "2025-01-30"
_TESTDB = "./testdb-integration"


@pytest.fixture(scope="session", autouse=True)
def setup_lamindb():
    ln.setup.init(storage=_TESTDB, modules="bionty")
    yield
    shutil.rmtree(_TESTDB, ignore_errors=True)
    ln.setup.delete("testdb-integration", force=True)


# ---------------------------------------------------------------------------
# Test 1: LTS smoke — registers artifacts with version_tag
# ---------------------------------------------------------------------------


def test_smoke_ingest_lts():
    """Smoke: registers the first 2 h5ad files from the LTS S3 path."""
    import cellxgene_lamin._register_annotate_new_release as _mod

    _mod.ingest_lts(new=LTS_NEW, previous=LTS_PREVIOUS, smoke=True)

    lts_artifacts = ln.Artifact.filter(version_tag=LTS_NEW)
    assert lts_artifacts.count() > 0, (
        f"expected LTS artifacts, got {lts_artifacts.count()}"
    )


# ---------------------------------------------------------------------------
# Test 2: pre-release smoke — LTS datasets skipped, new ones registered
# ---------------------------------------------------------------------------


def test_smoke_ingest_pre_release():
    """Smoke: LTS datasets are skipped; up to 2 non-LTS datasets are registered."""
    import cellxgene_lamin._register_annotate_new_release as _mod

    lts_dataset_ids = {
        af.key.split("/")[-1].replace(".h5ad", "")
        for af in ln.Artifact.filter(version_tag=LTS_NEW)
    }

    _mod.ingest_pre_release(new=LTS_NEW, smoke=True)

    pre_release_label = ln.ULabel.filter(name="pre-release").one_or_none()
    assert pre_release_label is not None, "pre-release ULabel was not created"

    pre_release_artifacts = ln.Artifact.filter(ulabels=pre_release_label)
    assert pre_release_artifacts.count() > 0, "no pre-release artifacts were registered"

    # none of the pre-release artifacts should be one of the LTS datasets
    for af in pre_release_artifacts:
        dataset_id = af.key.split("/")[-1].replace(".h5ad", "")
        assert dataset_id not in lts_dataset_ids, (
            f"dataset {dataset_id} is in LTS but was registered as pre-release"
        )


# ---------------------------------------------------------------------------
# Test 3: annotation links tissue and cell_type labels to the artifact
# ---------------------------------------------------------------------------


def test_annotation_links_tissue_and_cell_type_labels():
    """Annotation: after curating a pre-release artifact, tissue and cell_type labels are linked."""
    import cellxgene_lamin._register_annotate_new_release as _mod
    from cellxgene_lamin.dev._cxg_rest_api import get_datasets_from_cxg

    pre_release_label = ln.ULabel.filter(name="pre-release").one_or_none()
    assert pre_release_label is not None, "run test_smoke_ingest_pre_release first"

    afs = list(ln.Artifact.filter(ulabels=pre_release_label))
    assert afs, "no pre-release artifacts found to annotate"

    # annotate only the first artifact
    af = afs[0]
    dataset_id = af.key.split("/")[-1].replace(".h5ad", "")

    cxg_datasets = get_datasets_from_cxg()
    _mod._annotate_artifacts(
        cxg_datasets=cxg_datasets,
        registered_ids={dataset_id},
        new_census_version=LTS_NEW,
        pre_release_label=pre_release_label,
    )

    # re-fetch from DB after annotation
    af = ln.Artifact.get(uid=af.uid)
    assert af.cell_types.count() > 0, "no cell_type labels linked after annotation"
    assert af.tissues.count() > 0, "no tissue labels linked after annotation"
