from collections.abc import Iterable
from typing import Optional

from bionty.models import PublicSource, Registry

from ._features import FEATURE_TO_ACCESSOR, OBS_FEATURES


def create_ontology_record_from_source(
    ontology_id: str,
    from_orm: Registry,
    target_orm: Registry,
    public_source: Optional[PublicSource] = None,
):
    from_record = from_orm.from_public(
        ontology_id=ontology_id, public_source=public_source
    )
    try:
        target_record = target_orm(
            name=from_record.name,
            description=from_record.description,
            ontology_id=from_record.ontology_id,
            public_source_id=from_record.public_source_id,
        )
        return target_record
    except Exception:
        pass


def register_ontology_ids(cxg_datasets: Iterable):
    import bionty as bt
    import lamindb as ln

    ontology_ids = {}
    for name in OBS_FEATURES.keys():
        if name in ["donor_id", "suspension_type", "tissue_type"]:
            continue
        allids = set()
        for i in cxg_datasets:
            if name in i:
                allids.update([(j["label"], j["ontology_term_id"]) for j in i[name]])

        ontology_ids[name] = allids

    public_source_ds_mouse = bt.PublicSource.filter(
        entity="DevelopmentalStage", organism="mouse"
    ).one()
    public_source_pato = bt.PublicSource.filter(source="pato").one()

    upon_create_search_names = ln.settings.creation.search_names
    ln.settings.creation.search_names = False

    # register all ontology ids
    for name, terms in ontology_ids.items():
        print(f"registering {name}")
        accessor, orm = FEATURE_TO_ACCESSOR.get(name)
        terms_ids = [i[1] for i in terms]
        records = orm.from_values(terms_ids, field="ontology_id")
        if len(records) > 0:
            ln.save(records)
        inspect_result = orm.inspect(terms_ids, field="ontology_id", mute=True)
        if len(inspect_result.non_validated) > 0:
            if name == "development_stage":
                records = orm.from_values(
                    inspect_result.non_validated,
                    field="ontology_id",
                    public_source=public_source_ds_mouse,
                )
                records += [
                    create_ontology_record_from_source(
                        ontology_id=term_id, from_orm=bt.Tissue, target_orm=orm
                    )
                    for term_id in inspect_result.non_validated
                    if term_id.startswith("UBERON:")
                ]
                records += [
                    orm(name=term_id, ontology_id=term_id)
                    for term_id in inspect_result.non_validated
                    if term_id == "unknown"
                ]
            else:
                records = [
                    orm(name=term[0], ontology_id=term[1])
                    for term in terms
                    if (not term[1].startswith("PATO:"))
                    and (term[1] in inspect_result.non_validated)
                ]
                records += [
                    create_ontology_record_from_source(
                        ontology_id=term_id,
                        from_orm=bt.Phenotype,
                        target_orm=orm,
                        public_source=public_source_pato,
                    )
                    for term_id in inspect_result.non_validated
                    if term_id.startswith("PATO:")
                ]

            if len(records) > 0:
                print(f"registered {len(records)} records: {records}")
                ln.save(records)
    ln.settings.creation.search_names = upon_create_search_names

    # clean up the 2 "unknowns" in DevelopmentalStage
    bt.DevelopmentalStage.filter(name="unknown").exclude(ontology_id="unknown").delete()
