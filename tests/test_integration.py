"""Integration smoke tests — run against laminlabs/lamindata (real network calls)."""

import lamindb as ln
import pytest
from cellxgene_lamin._register_annotate_new_release import (
    ingest_lts,
    ingest_pre_release,
)

LTS_NEW = "2025-11-08"
LTS_PREVIOUS = "2025-01-30"


@pytest.fixture(scope="module", autouse=True)
def connect_lamindb():
    ln.connect("laminlabs/lamindata")


def test_smoke_ingest_pre_release():
    """Smoke: registers at most 2 pre-release artifacts against laminlabs/lamindata."""
    ingest_pre_release(new=LTS_NEW, smoke=True)

    pre_release_label = ln.ULabel.filter(name="pre-release").one_or_none()
    assert pre_release_label is not None

    registered = ln.Artifact.filter(ulabels=pre_release_label)
    assert registered.count() <= 2


def test_smoke_ingest_lts():
    """Smoke: registers at most 2 LTS artifacts against laminlabs/lamindata."""
    ingest_lts(new=LTS_NEW, previous=LTS_PREVIOUS, smoke=True)

    registered = ln.Artifact.filter(key__contains=LTS_NEW)
    assert registered.count() <= 2
