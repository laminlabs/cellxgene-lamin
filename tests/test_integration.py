"""Integration smoke tests — run against a temporary local lamindb instance."""

import shutil
import sys
from pathlib import Path

import lamindb as ln
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cellxgene_lamin._register_annotate_new_release import (
    ingest_lts,
    ingest_pre_release,
)

LTS_NEW = "2025-11-08"
LTS_PREVIOUS = "2025-01-30"


@pytest.fixture(scope="session", autouse=True)
def setup_lamindb():
    ln.setup.init(storage="./testdb")
    yield
    shutil.rmtree("./testdb")
    ln.setup.delete("testdb", force=True)


# ---------------------------------------------------------------------------
# Test 1: LTS smoke — registers artifacts with version_tag
# ---------------------------------------------------------------------------


def test_smoke_ingest_lts():
    """Smoke: registers the first 2 h5ad files from the LTS S3 path."""
    ingest_lts(new=LTS_NEW, previous=LTS_PREVIOUS, smoke=True)

    lts_artifacts = ln.Artifact.filter(version_tag=LTS_NEW)
    assert lts_artifacts.count() > 0, (
        f"expected LTS artifacts, got {lts_artifacts.count()}"
    )


# ---------------------------------------------------------------------------
# Test 2: pre-release smoke — LTS datasets skipped, new ones registered
# ---------------------------------------------------------------------------


def test_smoke_ingest_pre_release():
    """Smoke: LTS datasets are skipped; up to 2 non-LTS datasets are registered."""
    lts_dataset_ids = {
        af.key.split("/")[-1].replace(".h5ad", "")
        for af in ln.Artifact.filter(version_tag=LTS_NEW)
    }

    ingest_pre_release(new=LTS_NEW, smoke=True)

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
