import pandas as pd
from lnschema_core.types import FieldAttr


def convert_name_to_ontology_id(values: pd.Series, field: FieldAttr):
    field_name = field.field.name
    assert field_name == "name"
    cols = ["name", "ontology_id"]
    registry = field.field.model
    if hasattr(registry, "ontology_id"):
        validated_records = registry.filter(**{f"{field_name}__in": values})
        mapper = (
            pd.DataFrame(validated_records.values_list(*cols)).set_index(0).to_dict()[1]
        )
        return values.map(mapper)
