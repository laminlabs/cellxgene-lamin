import shutil
from pathlib import Path

import nox
from laminci import convert_executable_md_files, upload_docs_artifact
from laminci.nox import (
    build_docs,
    install_lamindb,
    login_testuser1,
    run,
    run_pre_commit,
)

nox.options.default_venv_backend = "none"

GROUPS = {
    "curate": ["cellxgene-curate.ipynb"],
    "query": ["cellxgene.ipynb"],
}


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session
@nox.parametrize(
    "group",
    ["query", "curate", "docs"],
)
def install(session: nox.Session, group: str) -> None:
    extras = ""
    if group != "docs":
        run(
            session,
            "uv pip install --system pronto tiledbsoma scanpy>=1.11.3",
        )  # scanpy pin to prevent scipy installation crashes
        run(session, "uv tool install cellxgene-schema==5.2.2")
    install_lamindb(session, branch="main", extras=extras)
    run(session, "uv pip install --system .[dev]")


@nox.session
@nox.parametrize(
    "group",
    ["query", "curate"],
)
def build(session, group):
    convert_executable_md_files()
    login_testuser1(session)
    if group != "curate":
        run(session, f"pytest -s ./tests/test_notebooks.py::test_{group}")

    # Move executed notebooks temporarily to recover them for docs building
    target_dir = Path(f"./docs_{group}")
    target_dir.mkdir(exist_ok=True)
    for filename in GROUPS[group]:
        shutil.copy(Path("docs") / filename, target_dir / filename)


@nox.session
def docs(session):
    convert_executable_md_files()
    # Recover executed notebooks
    for group in ["query", "curate"]:
        for path in Path(f"./docs_{group}").glob("*"):
            path.rename(f"./docs/{path.name}")

    run(session, "lamin init --storage ./docsbuild --modules bionty")
    build_docs(session, strict=False)
    upload_docs_artifact()
