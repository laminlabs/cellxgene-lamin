import sys
from pathlib import Path

import nbproject_test as test

sys.path[:0] = [str(Path(__file__).parent.parent)]

from noxfile import GROUPS  # noqa: E402

DOCS = Path(__file__).parents[1] / "docs/"


def test_census():
    for filename in GROUPS["census"]:
        print(filename)
        test.execute_notebooks(DOCS / filename, write=True)


def test_validator():
    for filename in GROUPS["validator"]:
        print(filename)
        test.execute_notebooks(DOCS / filename, write=True)
