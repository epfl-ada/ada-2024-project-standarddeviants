import json
import os
from typing import Iterable

import pandas as pd
import requests
import tqdm


def get_uniprot_diseases(
    uniprot_ids: Iterable, savename: str | os.PathLike = None
) -> dict[str, dict[str, list[str]]]:
    """Get the diseases mentionned in the Uniprot Database from the ID of a target chain.

    Args:
        uniprot_ids (Iterable): Target chain IDs used to find the related diseases.
        savename (str | os.PathLike, optional): The path to save the output to in json format. Defaults to None.

    Returns:
        dict[str, dict[str, list[str]]]: A dictionaary of the form : { Uniprot_ID: { "comments": [...],  "keywords": [...] }, ... } where each list enumerates the related diseases.
    """
    diseases = {}
    for id in tqdm.tqdm(uniprot_ids):
        r = requests.get(f"https://www.uniprot.org/uniprotkb/{id}.json").json()
        diseases_per_id = {
            "comments": [],
            "keywords": [],
        }
        for comment in r.get("comments", []):
            if (comment["commentType"] == "DISEASE") and (
                comment.get("disease") is not None
            ):
                diseaseID = comment.get("disease").get("diseaseId")
                if diseaseID is not None:
                    diseases_per_id["comments"].append(diseaseID)

        for keyword in r.get("keywords", []):
            if keyword["category"].lower() == "disease":
                diseases_per_id["keywords"].append(keyword["name"])

        diseases.update({id: diseases_per_id})

    if savename is not None:
        with open(savename, "w") as f:
            json.dump(diseases, f)

    return diseases


if __name__ == "__main__":
    df = pd.read_csv(
        r"../data/BindingDB_All.tsv",
        sep="\t",
        usecols=["UniProt (SwissProt) Primary ID of Target Chain"],
    )

    get_uniprot_diseases(
        df["UniProt (SwissProt) Primary ID of Target Chain"].value_counts().index,
        savename=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "data",
            "UniprotID_disases.json",
        ),
    )
