import lamindb as ln
from lamindb.core._feature_manager import get_accessor_by_registry_

OBS_FEATURES = {
    "assay": "bionty.ExperimentalFactor",
    "cell_type": "bionty.CellType",
    "development_stage": "bionty.DevelopmentalStage",
    "disease": "bionty.Disease",
    "donor_id": "core.ULabel",
    "self_reported_ethnicity": "bionty.Ethnicity",
    "sex": "bionty.Phenotype",
    "suspension_type": "core.ULabel",
    "tissue": "bionty.Tissue",
    "tissue_type": "core.ULabel",
    "organism": "bionty.Organism",
}

obs_features_records = ln.FeatureSet.filter(name="obs metadata").one().members.lookup()
ACCESSORS = get_accessor_by_registry_(ln.Artifact)
FEATURE_TO_ACCESSOR = {}
for name in OBS_FEATURES.keys():
    feature = getattr(obs_features_records, name)
    accessor = ACCESSORS.get(
        feature.dtype[feature.dtype.index("[") + 1 : feature.dtype.index("]")]
    )
    orm = getattr(ln.Artifact, accessor).field.model
    # TODO: ulabels are defined in the File model, improve this in LaminDB
    if orm == ln.Artifact:
        orm = getattr(ln.Artifact, accessor).field.related_model
    FEATURE_TO_ACCESSOR[name] = (accessor, orm)


def register_obs_featureset(artifacts):
    import lamindb as ln

    feature_set = ln.FeatureSet.filter(name="obs metadata").one_or_none()
    if feature_set is None:
        features_records = []
        for name, registry in OBS_FEATURES.items():
            record = ln.Feature(name=name, dtype=f"cat[{registry}]").save()
            features_records.append(record)
        ln.save(features_records)
        feature_set = ln.FeatureSet(features=features_records, name="obs metadata")
        feature_set.save()

    feature_set.artifacts.add(*artifacts, through_defaults={"slot": "obs"})
    return feature_set
