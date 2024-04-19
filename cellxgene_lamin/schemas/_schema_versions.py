from pathlib import Path

import pandas as pd
import yaml


def read_schema_versions(ontology_versions: Path) -> dict[str, pd.DataFrame]:
    data = yaml.safe_load(open(ontology_versions))
    schema_versions = data["schema-version"]

    def _schema_to_df(schema_data):
        rows = []
        for entity, details in schema_data.items():
            for ontology, values in details.items():
                for organism, version in values.items():
                    rows.append((entity, organism, ontology, version))
        return pd.DataFrame(
            rows, columns=["entity", "organism", "source", "version"]
        ).set_index("entity")

    desired_dfs_from_yaml = {
        version: _schema_to_df(details) for version, details in schema_versions.items()
    }

    return desired_dfs_from_yaml
