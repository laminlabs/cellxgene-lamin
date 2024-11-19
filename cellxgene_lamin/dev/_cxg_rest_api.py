from typing import Any

import requests


def get_datasets_from_cxg() -> dict[str, Any]:
    api_url_base = "https://api.cellxgene.cziscience.com"
    datasets_path = "/curation/v1/datasets"
    datasets_url = f"{api_url_base}{datasets_path}"
    headers = {"Content-Type": "application/json"}
    res = requests.get(url=datasets_url, headers=headers)
    res.raise_for_status()
    res_content = res.json()
    return res_content


def get_collections_from_cxg() -> dict[str, Any]:
    api_url_base = "https://api.cellxgene.cziscience.com"
    collections_path = "/curation/v1/collections"
    collections_url = f"{api_url_base}{collections_path}"
    headers = {"Content-Type": "application/json"}
    res = requests.get(url=collections_url, headers=headers)
    res.raise_for_status()
    res_content = res.json()
    return res_content
