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
        session.run(*"uv pip install --system cellxgene-census".split())
    elif group == "validator":
        extra = ",jupyter,aws,zarr"
        session.run(*"uv pip install --system cellxgene-schema".split())
        session.run(*"uv pip install --system anndata==0.9.0".split())
    session.run(*"uv pip install --system .[dev]".split())
    session.run(
        "pip",
        "install",
        f"lamindb[bionty{extra}] @ git+https://github.com/laminlabs/lamindb@main",
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


@nox.session
def docs(session):
    login_testuser1(session)
    extra = ",jupyter,aws,zarr"
    session.run(*"uv pip install --system cellxgene-schema".split())
    session.run(*"uv pip install --system anndata==0.9.0".split())
    session.run(*"uv pip install --system .[dev]".split())
    session.run(
        "pip",
        "install",
        f"lamindb[bionty{extra}] @ git+https://github.com/laminlabs/lamindb@main",
    )
    session.run(*"pytest -s ./tests/test_notebooks.py::test_validator".split())
    session.run(*"lamin init --storage ./docsbuild --schema bionty".split())
    build_docs(session, strict=True)
    upload_docs_artifact(aws=True)
