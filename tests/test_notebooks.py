from pathlib import Path

import nbproject_test as test
from noxfile import GROUPS

DOCS = Path(__file__).parent.parent / "docs/"


def test_census():
    for filename in GROUPS["census"]:
        print(filename)
        test.execute_notebooks(DOCS / filename, write=True)


def test_validator():
    for filename in GROUPS["validator"]:
        print(filename)
        test.execute_notebooks(DOCS / filename, write=True)
