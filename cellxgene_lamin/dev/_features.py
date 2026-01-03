import lamindb as ln
from lamindb.models._feature_manager import get_accessor_by_registry_

OBS_FEATURES = {
    "organism": "bionty.Organism",
    "assay": "bionty.ExperimentalFactor",
    "cell_type": "bionty.CellType",
    "development_stage": "bionty.DevelopmentalStage",
    "disease": "bionty.Disease",
    "self_reported_ethnicity": "bionty.Ethnicity",
    "sex": "bionty.Phenotype",
    "tissue": "bionty.Tissue",
    "suspension_type": "ULabel",
    "donor_id": "ULabel",
    "tissue_type": "ULabel",
}


class _FeatureAccessorMapping:
    def __init__(self):
        self._mapping = None

    def _compute_mapping(self):
        obs_features_records = (
            ln.Schema.filter(name="obs metadata").one().members.lookup()
        )
        accessors = get_accessor_by_registry_(ln.Artifact)
        mapping = {}

        for name in OBS_FEATURES.keys():
            feature = getattr(obs_features_records, name)
            accessor = accessors.get(
                feature._dtype_str[
                    feature._dtype_str.index("[") + 1 : feature._dtype_str.index("]")
                ]
            )
            orm = getattr(ln.Artifact, accessor).field.model
            if orm == ln.Artifact:
                orm = getattr(ln.Artifact, accessor).field.related_model
            mapping[name] = (accessor, orm)

        return mapping

    def __getitem__(self, key):
        if self._mapping is None:
            self._mapping = self._compute_mapping()
        return self._mapping[key]

    def __iter__(self):
        if self._mapping is None:
            self._mapping = self._compute_mapping()
        return iter(self._mapping)

    def get(self, key, default=None):
        if self._mapping is None:
            self._mapping = self._compute_mapping()
        return self._mapping.get(key, default)


FEATURE_TO_ACCESSOR = _FeatureAccessorMapping()


def register_obs_schema(artifacts: ln.Artifact) -> ln.Schema:
    import lamindb as ln

    feature_set = ln.Schema.filter(name="obs metadata").one_or_none()
    if feature_set is None:
        features_records = []
        for name, registry in OBS_FEATURES.items():
            record = ln.Feature(name=name, dtype=f"cat[{registry}]").save()
            features_records.append(record)
        ln.save(features_records)
        feature_set = ln.Schema(features=features_records, name="obs metadata")
        feature_set.save()

    feature_set.artifacts.add(*artifacts, through_defaults={"slot": "obs"})

    return feature_set
