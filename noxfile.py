import shutil
from pathlib import Path

import nox
from laminci import upload_docs_artifact
from laminci.nox import build_docs, login_testuser1, run_pre_commit

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
    ["census", "validator", "docs"],
)
def install(session: nox.Session, group: str) -> None:
    extra = ""
    session.run(*"uv pip install --system scipy>=1.12.0,<1.13.0rc1".split())
    if group == "census":
        extra = ",jupyter,aws"
        session.run(*"uv pip install --system tiledbsoma".split())
    elif group == "validator":
        extra = ",jupyter,aws,zarr"
        session.run(*"uv pip install --system cellxgene-schema==5.0.2".split())
    session.run(*"uv pip install --system .[dev]".split())
    session.run(
        "uv",
        "pip",
        "install",
        "--system",
        f"lamindb[bionty{extra}] @ git+https://github.com/laminlabs/lamindb@main",
    )
    session.run(
        "uv",
        "pip",
        "install",
        "--system",
        "git+https://github.com/laminlabs/bionty@main",
    )


@nox.session
@nox.parametrize(
    "group",
    ["census", "validator"],
)
@nox.session
def build(session, group):
    login_testuser1(session)
    session.run(*f"pytest -s ./tests/test_notebooks.py::test_{group}".split())

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

    session.run(*"lamin init --storage ./docsbuild --schema bionty".split())
    build_docs(session, strict=True)
    upload_docs_artifact(aws=True)
