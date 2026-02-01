from pathlib import Path

import nbproject_test as test

GROUPS = {
    "query": ["cellxgene.ipynb"],
    "curate": ["cellxgene-curate.ipynb"],
}
DOCS = Path(__file__).parent.parent / "docs/"


def test_query():
    for filename in GROUPS["query"]:
        test.execute_notebooks(DOCS / filename, write=True)


def test_validator():
    for filename in GROUPS["curate"]:
        test.execute_notebooks(DOCS / filename, write=True)
