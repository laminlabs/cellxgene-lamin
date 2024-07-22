from pathlib import Path

import nbproject_test as test

GROUPS = {
    "census": ["query-census.ipynb"],
    "validator": ["cellxgene.ipynb", "cellxgene-curate.ipynb"],
}
DOCS = Path(__file__).parent.parent / "docs/"


def test_census():
    for filename in GROUPS["census"]:
        test.execute_notebooks(DOCS / filename, write=True)


def test_validator():
    for filename in GROUPS["validator"]:
        test.execute_notebooks(DOCS / filename, write=True)
