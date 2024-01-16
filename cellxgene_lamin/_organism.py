def register_organisms(cxg_datasets):
    import lamindb as ln
    import lnschema_bionty as lb

    # register all organisms
    ncbitaxon_source = lb.BiontySource.filter(source="ncbitaxon").one()

    organisms_meta = set()
    for cxg_dataset in cxg_datasets:
        organisms_meta.update({i["ontology_term_id"] for i in cxg_dataset["organism"]})

    organisms_records = lb.Organism.from_values(
        organisms_meta, field=lb.Organism.ontology_id, bionty_source=ncbitaxon_source
    )
    # rename house mouse to mouse
    for r in organisms_records:
        if r.name == "house mouse":
            r.name = "mouse"
    ln.save(organisms_records, parents=False)


def annotate_organisms(artifacts, cxg_datasets):
    import lamindb as ln
    import lnschema_bionty as lb

    ext_feature_set = ln.FeatureSet.filter(name="external metadata").one()
    ext_features = ext_feature_set.members.lookup()

    for cxg_dataset in cxg_datasets:
        # get registered file record based on dataset_id
        artifact = artifacts.filter(
            key__contains=cxg_dataset["dataset_id"]
        ).one_or_none()
        if artifact is None:
            continue

        # annotate artifacts with organisms
        organism_ontology_ids = [i["ontology_term_id"] for i in cxg_dataset["organism"]]
        organism_records = lb.Organism.filter(
            ontology_id__in=organism_ontology_ids
        ).list()
        artifact.labels.add(organism_records, feature=ext_features.organism)
