import json
import time
from collections.abc import Iterable

import requests


def get_ensembl_versions(ensembl_ids: Iterable[str]) -> dict[str, int | None]:
    base_url = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json", "User-Agent": "anonymous"}
    results = {}

    for ensembl_id in ensembl_ids:
        endpoint = f"/archive/id/{ensembl_id}"
        try:
            response = requests.get(base_url + endpoint, headers=headers)
            response.raise_for_status()
            data = response.json()
            release = data.get("release")
            if release:
                results[ensembl_id] = int(release)
            else:
                print(f"No release found for {ensembl_id}")
                results[ensembl_id] = None
        except requests.exceptions.RequestException as e:
            print(f"Error for {ensembl_id}: {str(e)}")
            if hasattr(e, "response"):
                print(f"Error response: {e.response.text}")
            results[ensembl_id] = None
        except json.JSONDecodeError:
            print(f"Error decoding JSON for {ensembl_id}")
            print(f"Response content: {response.text}")
            results[ensembl_id] = None

        time.sleep(0.1)  # To ensure that we don't run into rate limiting

    return results


ensembl_ids = ["ENSG00000139618", "ENSG00000141510", "ENSG00000157764"]
results = get_ensembl_versions(ensembl_ids)

for ensembl_id, release in results.items():
    if release:
        print(f"Ensembl ID: {ensembl_id}, Release: {release}")
    else:
        print(f"Ensembl ID: {ensembl_id}, Unable to fetch release")
