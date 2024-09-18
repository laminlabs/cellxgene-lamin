from __future__ import annotations

import re
from importlib import resources
from typing import TYPE_CHECKING, Literal

import bionty as bt
import pandas as pd
from lamin_utils import logger
from lamindb._curate import AnnDataCurator, validate_categories_in_df

from .fields import CellxGeneFields
from .schemas._schema_versions import _read_schema_versions

if TYPE_CHECKING:
    from pathlib import Path

    import anndata as ad
    from anndata import AnnData
    from lnschema_core import Record
    from lnschema_core.types import FieldAttr


def _restrict_obs_fields(
    adata: ad.AnnData, obs_fields: dict[str, FieldAttr]
) -> dict[str, str]:
    """Restrict the obs fields to name return only available obs fields.

    To simplify the curation, we only validate against either name or ontology_id.
    If both are available, we validate against ontology_id.
    If none are available, we validate against name.
    """
    obs_fields_unique = {k: v for k, v in obs_fields.items() if k in adata.obs.columns}
    for name, field in obs_fields.items():
        if name.endswith("_ontology_term_id"):
            continue
        # if both the ontology id and the name are present, only validate on the ontology_id
        if (
            name in adata.obs.columns
            and f"{name}_ontology_term_id" in adata.obs.columns
        ):
            obs_fields_unique.pop(name)
        # if the neither name nor ontology id are present, validate on the name
        # this will raise error downstream, we just use name to be more readable
        if (
            name not in adata.obs.columns
            and f"{name}_ontology_term_id" not in adata.obs.columns
        ):
            obs_fields_unique[name] = field

    # Only retain obs_fields_unique that have keys in adata.obs.columns
    available_obs_fields = {
        k: v for k, v in obs_fields_unique.items() if k in adata.obs.columns
    }

    return available_obs_fields


def _add_defaults_to_obs(
    adata: ad.AnnData,
    defaults: dict[str, str],
) -> None:
    """Add defaults to the obs fields."""
    added_defaults: dict = {}
    for name, default in defaults.items():
        if (
            name not in adata.obs.columns
            and f"{name}_ontology_term_id" not in adata.obs.columns
        ):
            adata.obs[name] = default
            added_defaults[name] = default
    if len(added_defaults) > 0:
        logger.important(f"added defaults to the AnnData object: {added_defaults}")


class Curate(AnnDataCurator):
    """Annotation flow of AnnData based on CELLxGENE schema."""

    def __init__(
        self,
        adata: ad.AnnData | str | Path,
        var_index: FieldAttr = bt.Gene.ensembl_gene_id,
        categoricals: dict[str, FieldAttr] = CellxGeneFields.OBS_FIELDS,
        *,
        defaults: dict[str, str] = None,
        using_key: str = "laminlabs/cellxgene",
        verbosity: str = "hint",
        organism: str | None = None,
        schema_version: Literal["4.0.0", "5.0.0", "5.1.0"] = "5.1.0",
    ):
        self.organism = organism
        self.using_key = using_key

        VALID_SCHEMA_VERSIONS = {"4.0.0", "5.0.0", "5.1.0"}
        if schema_version not in VALID_SCHEMA_VERSIONS:
            valid_versions = ", ".join(sorted(VALID_SCHEMA_VERSIONS))
            raise ValueError(
                f"Invalid schema_version: {schema_version}. "
                f"Valid versions are: {valid_versions}"
            )
        self.schema_version = schema_version
        self.schema_reference = f"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/{schema_version}/schema.md"
        with resources.path(
            "cellxgene_lamin.schemas", "schema_versions.yml"
        ) as schema_versions_path:
            self._pinned_ontologies = _read_schema_versions(schema_versions_path)[
                self.schema_version
            ]
        self.sources = self._create_sources(adata)
        self.sources = {
            entity: source
            for entity, source in self.sources.items()
            if source is not None
        }

        if defaults:
            _add_defaults_to_obs(adata, defaults)

        exclude_keys = {
            entity: default
            for entity, default in CellxGeneFields.OBS_FIELD_DEFAULTS.items()
            if entity in adata.obs.columns  # type: ignore
        }
        super().__init__(
            data=adata,
            var_index=var_index,
            categoricals=_restrict_obs_fields(adata, categoricals),
            using_key=using_key,
            verbosity=verbosity,
            organism=organism,
            sources=self.sources,
            exclude=exclude_keys,
        )

    @property
    def pinned_ontologies(self) -> pd.DataFrame:
        print(f"Currently used schema version: {self.schema_version}")
        return self._pinned_ontologies

    @property
    def adata(self) -> AnnData:
        return self._adata

    def _create_sources(self, adata: ad.AnnData) -> dict[str, Record]:
        """Creates a sources dictionary that can be passed to AnnDataCurator."""

        # fmt: off
        def _fetch_bionty_source(
            entity: str, organism: str, source: str
        ) -> bt.Source | None:
            """Fetch the Bionty source of the pinned ontology.

            Returns None if the source does not exist.
            """
            version = self._pinned_ontologies.loc[(self._pinned_ontologies.index == entity) &
                                                  (self._pinned_ontologies["organism"] == organism) &
                                                  (self._pinned_ontologies["source"] == source), "version"].iloc[0]
            return bt.Source.using(self.using_key).filter(organism=organism, entity=f"bionty.{entity}", version=version).first()

        entity_mapping = {
             "var_index": ("Gene", self.organism, "ensembl"),
             "gene": ("Gene", self.organism, "ensembl"),
             "cell_type": ("CellType", "all", "cl"),
             "assay": ("ExperimentalFactor", "all", "efo"),
             "self_reported_ethnicity": ("Ethnicity", self.organism, "hancestro"),
             "development_stage": ("DevelopmentalStage", self.organism, "hsapdv" if self.organism == "human" else "mmusdv"),
             "disease": ("Disease", "all", "mondo"),
             # "organism": ("Organism", "vertebrates", "ensembl"),
             "sex": ("Phenotype", "all", "pato"),
             "tissue": ("Tissue", "all", "uberon"),
        }
        # fmt: on

        entity_to_sources = {
            key: source
            for entity, params in entity_mapping.items()
            for key, source in {
                entity: _fetch_bionty_source(*params),
                f"{entity}_ontology_id": _fetch_bionty_source(*params),
            }.items()
            if key in adata.obs.columns or key == "var_index"
        }

        return entity_to_sources

    def _convert_name_to_ontology_id(self, values: pd.Series, field: FieldAttr):
        """Converts a column that stores a name into a column that stores the ontology id.

        cellxgene expects the obs columns to be {entity}_ontology_id columns and disallows {entity} columns.
        """
        field_name = field.field.name
        assert field_name == "name"
        cols = ["name", "ontology_id"]
        registry = field.field.model

        if hasattr(registry, "ontology_id"):
            validated_records = registry.using(self.using_key).filter(
                **{f"{field_name}__in": values}
            )
            mapper = (
                pd.DataFrame(validated_records.values_list(*cols))
                .set_index(0)
                .to_dict()[1]
            )
            return values.map(mapper)

    def validate(self) -> bool:
        """Validates the AnnData object against most cellxgene requirements."""
        # Verify that all required obs columns are present
        missing_obs_fields = [
            name
            for name in CellxGeneFields.OBS_FIELD_DEFAULTS.keys()
            if name not in self._adata.obs.columns
            and f"{name}_ontology_term_id" not in self._adata.obs.columns
        ]

        if len(missing_obs_fields) > 0:
            missing_obs_fields_str = ", ".join(list(missing_obs_fields))
            logger.error(f"missing required obs columns {missing_obs_fields_str}")
            logger.info(
                "consider initializing a Curate object like 'Curate(adata, defaults=cxg.CellxGeneFields.OBS_FIELD_DEFAULTS)'"
                "to automatically add these columns with default values."
            )
            return False

        # Verify that no cellxgene reserved names are present
        reserved_names = {
            "ethnicity",
            "ethnicity_ontology_term_id",
            "X_normalization",
            "default_field",
            "layer_descriptions",
            "tags",
            "versions",
            "contributors",
            "preprint_doi",
            "project_description",
            "project_links",
            "project_name",
            "publication_doi",
        }
        matched_columns = [
            column for column in self._adata.obs.columns if column in reserved_names
        ]
        if len(matched_columns) > 0:
            raise ValueError(
                f"AnnData object must not contain obs columns {matched_columns} which are"
                " reserved from previous schema versions."
            )

        # cellxgene requires an embedding
        embedding_pattern = r"^[a-zA-Z][a-zA-Z0-9_.-]*$"
        exclude_key = "spatial"
        matching_keys = [
            key
            for key in self._adata.obsm.keys()
            if re.match(embedding_pattern, key) and key != exclude_key
        ]
        if len(matching_keys) == 0:
            raise ValueError(
                "Unable to find an embedding key. Please calculate an embedding."
            )

        return super().validate(organism=self.organism)

    def to_cellxgene_anndata(
        self, is_primary_data: bool, title: str | None = None
    ) -> ad.AnnData:
        """Converts the AnnData object to the cellxgene-schema input format.

        cellxgene expects the obs fields to be {entity}_ontology_id fields and has many further requirements which are
        documented here: https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema.
        This function checks for most but not all requirements of the CELLxGENE schema.
        If you want to ensure that it fully adheres to the CELLxGENE schema, run `cellxgene-schema` on the AnnData object.

        Args:
            is_primary_data: Whether the measured data is primary data or not.
            title: Title of the AnnData object. Commonly the name of the publication.

        Returns:
            An AnnData object which adheres to the cellxgene-schema.
        """
        # Create a copy since we modify the AnnData object extensively
        adata_cxg = self._adata.copy()

        # convert name column to ontology_term_id column
        for column in adata_cxg.obs.columns:
            if column in self.categoricals and not column.endswith("_ontology_term_id"):
                mapped_column = self._convert_name_to_ontology_id(
                    adata_cxg.obs[column], field=self.categoricals.get(column)
                )
                if mapped_column is not None:
                    adata_cxg.obs[f"{column}_ontology_term_id"] = mapped_column

        # drop the name columns for ontologies. cellxgene does not allow them.
        drop_columns = [
            i
            for i in adata_cxg.obs.columns
            if f"{i}_ontology_term_id" in adata_cxg.obs.columns
        ]
        adata_cxg.obs.drop(columns=drop_columns, inplace=True)

        # Add cellxgene metadata to AnnData object
        if "is_primary_data" not in adata_cxg.obs.columns:
            adata_cxg.obs["is_primary_data"] = is_primary_data
        if "feature_is_filtered" not in adata_cxg.var.columns:
            logger.warn(
                "column 'feature_is_filtered' not present in var. Setting to default"
                " value of False."
            )
            adata_cxg.var["feature_is_filtered"] = False
        if self._collection is None:
            if title is None:
                raise ValueError("please pass a title!")
            else:
                adata_cxg.uns["title"] = title
        else:
            adata_cxg.uns["title"] = self._collection.name
        adata_cxg.uns["cxg_lamin_schema_reference"] = self.schema_reference
        adata_cxg.uns["cxg_lamin_schema_version"] = self.schema_version

        return adata_cxg
