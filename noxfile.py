import shutil
from pathlib import Path

import nox
from laminci import upload_docs_artifact
from laminci.nox import (
    build_docs,
    install_lamindb,
    login_testuser1,
    run,
    run_pre_commit,
)

nox.options.default_venv_backend = "none"

GROUPS = {
    "census": ["query-census.ipynb"],
    "validator": ["cellxgene.ipynb", "cellxgene-curate.ipynb"],
}


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session
@nox.parametrize(
    "group",
    ["validator", "docs"],
)
def install(session: nox.Session, group: str) -> None:
    extras = ""
    if group == "validator":
        extras = "bionty,jupyter,zarr"
        run(session, "uv pip install --system tiledbsoma==1.15.0rc3")
        run(session, "uv tool install cellxgene-schema==5.2.2")
    install_lamindb(session, branch="main", extras=extras)
    run(session, "uv pip install --system .[dev]")


@nox.session
@nox.parametrize(
    "group",
    ["census", "validator"],
)
@nox.session
def build(session, group):
    login_testuser1(session)
    run(session, f"pytest -s ./tests/test_notebooks.py::test_{group}")

    # Move executed notebooks temporarily to recover them for docs building
    target_dir = Path(f"./docs_{group}")
    target_dir.mkdir(exist_ok=True)
    for filename in GROUPS[group]:
        shutil.copy(Path("docs") / filename, target_dir / filename)


@nox.session
def docs(session):
    # Recover executed notebooks
    for group in ["census", "validator"]:
        for path in Path(f"./docs_{group}").glob("*"):
            path.rename(f"./docs/{path.name}")

    run(session, "lamin init --storage ./docsbuild --schema bionty")
    build_docs(session, strict=True)
    upload_docs_artifact(aws=True)
