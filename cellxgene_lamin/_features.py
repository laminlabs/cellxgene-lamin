import lamindb as ln
from lamindb.dev._feature_manager import get_accessor_by_orm

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
}

EXT_FEATURES = {"organism": "bionty.Organism"}

obs_features_records = ln.FeatureSet.filter(name="obs metadata").one().members.lookup()
ACCESSORS = get_accessor_by_orm(ln.Artifact)
FEATURE_TO_ACCESSOR = {}
for name in OBS_FEATURES.keys():
    feature = getattr(obs_features_records, name)
    accessor = ACCESSORS.get(feature.registries)
    orm = getattr(ln.Artifact, accessor).field.model
    # TODO: ulabels are defined in the File model, improve this in LaminDB
    if orm == ln.Artifact:
        orm = getattr(ln.Artifact, accessor).field.related_model
    FEATURE_TO_ACCESSOR[name] = (accessor, orm)


def register_feature_set(artifacts, slot: str):
    import lamindb as ln

    if slot == "obs":
        features = OBS_FEATURES
    elif slot == "ext":
        features = EXT_FEATURES

    feature_set = ln.FeatureSet.filter(name=f"{slot} features").one_or_none()
    if feature_set is None:
        features_records = []
        for name, registry in features.items():
            record = ln.Feature(name=name, type="category", registries=registry)
            features_records.append(record)
        ln.save(features_records)
        feature_set = ln.FeatureSet(features=features_records, name=f"{slot} metadata")
        feature_set.save()

    feature_set.artifacts.add(*artifacts, through_defaults={"slot": slot})
    return feature_set
