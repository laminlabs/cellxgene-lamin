import shutil
import sys
from pathlib import Path

import lamindb as ln
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

_TESTDB = "./testdb-integration"


@pytest.fixture(scope="session", autouse=True)
def setup_lamindb():
    ln.setup.init(storage=_TESTDB, modules="bionty")
    yield
    shutil.rmtree(_TESTDB, ignore_errors=True)
    ln.setup.delete("testdb-integration", force=True)


@pytest.fixture(scope="session")
def mod(setup_lamindb):
    # imported here — after ln.setup.init — so Django is fully initialized
    # before lamindb models in the module are loaded
    import cellxgene_lamin._register_annotate_new_release as _mod

    return _mod
