import bionty as bt
import lamindb as ln


class CellxGeneFields:
    """CELLxGENE fields."""

    OBS_FIELDS = {
        "assay": bt.ExperimentalFactor.name,
        "assay_ontology_term_id": bt.ExperimentalFactor.ontology_id,
        "cell_type": bt.CellType.name,
        "cell_type_ontology_term_id": bt.CellType.ontology_id,
        "development_stage": bt.DevelopmentalStage.name,
        "development_stage_ontology_term_id": bt.DevelopmentalStage.ontology_id,
        "disease": bt.Disease.name,
        "disease_ontology_term_id": bt.Disease.ontology_id,
        "donor_id": ln.ULabel.name,
        "self_reported_ethnicity": bt.Ethnicity.name,
        "self_reported_ethnicity_ontology_term_id": bt.Ethnicity.ontology_id,
        "sex": bt.Phenotype.name,
        "sex_ontology_term_id": bt.Phenotype.ontology_id,
        "suspension_type": ln.ULabel.name,
        "tissue": bt.Tissue.name,
        "tissue_ontology_term_id": bt.Tissue.ontology_id,
        "tissue_type": ln.ULabel.name,
        "organism": bt.Organism.name,
        "organism_ontology_term_id": bt.Organism.ontology_id,
    }

    OBS_FIELD_DEFAULTS = {
        "disease": "normal",
        "development_stage": "unknown",
        "self_reported_ethnicity": "unknown",
        "suspension_type": "cell",
        "donor_id": "na",
        "tissue_type": "tissue",
        "cell_type": "native_cell",
        "sex": "unknown",
    }
