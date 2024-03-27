from pathlib import Path

import nbproject_test

notebook_dir = Path(__file__).parent.parent / "docs/"


def test_all_notebooks():
    nbproject_test.execute_notebooks(notebook_dir)
