from collections.abc import Iterable


def register_ulabels(cxg_datasets: Iterable, feature_name: str):
    import lamindb as ln

    if feature_name in {"donor_id", "suspension_type"}:
        labels = set()
        for i in cxg_datasets:
            if feature_name in i:
                labels.update(i[feature_name])

        is_feature_name = ln.ULabel.filter(name=f"is_{feature_name}").one_or_none()
        if is_feature_name is None:
            is_feature_name = ln.ULabel(
                name=f"is_{feature_name}", description=f"parent of {feature_name}s"
            ).save()
        records = is_feature_name.children.all()
        result = records.inspect(labels, mute=True)
        new_records = [ln.ULabel(name=name) for name in result.non_validated]
        ln.save(new_records)
        print(f"registered {len(new_records)} {feature_name}s")
        is_feature_name.children.add(*new_records)
    elif feature_name == "tissue_type":
        is_tissue_type = ln.ULabel.filter(name="is_tissue_type").one_or_none()
        if is_tissue_type is None:
            is_tissue_type = ln.ULabel(
                name="is_tissue_type", description="parents of tissue types"
            )
            is_tissue_type.save()
            tissue_types = [
                ln.ULabel(name=i) for i in ["tissue", "organoid", "cell culture"]
            ]
            ln.save(tissue_types)
            is_tissue_type.children.add(*tissue_types)
