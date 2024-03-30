import shutil
from pathlib import Path

import nox
from laminci import upload_docs_artifact
from laminci.nox import build_docs, login_testuser1, run_pre_commit

nox.options.default_venv_backend = "none"

GROUPS = {
    "census": ["query-census.ipynb"],
    "validator": ["cellxgene.ipynb", "cellxgene-lamin-validator.ipynb"],
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
    if group == "census":
        extra = ",jupyter,aws"
        session.run(*"pip install cellxgene-census".split())
    elif group == "validator":
        extra = ",jupyter,aws,zarr"
        session.run(*"pip install cellxgene-schema".split())
        session.run(*"pip install anndata==0.9.0".split())
    session.run(*"pip install .[dev]".split())
    session.run(
        "pip",
        "install",
        f"lamindb[bionty{extra}] @ git+https://github.com/laminlabs/lamindb@release",
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

    login_testuser1(session)
    session.run(*"lamin init --storage ./docsbuild --schema bionty".split())
    build_docs(session, strict=False)
    upload_docs_artifact(aws=True)
