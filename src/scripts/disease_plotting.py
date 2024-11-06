import json
import os

import numpy as np
import pandas as pd


def load_uniprotid_diseases(
    filepath: str | os.PathLike,
) -> dict[str, dict[str, list[str]]]:
    """Load the dictionnary saved by the `get_uniprot_diseases` function"""
    with open(filepath, "r") as f:
        diseases = json.load(f)
    return diseases


def quantify_missing_diseases(
    diseases: dict[str, dict[str, list[str]]]
) -> tuple[list[str], float]:
    """Returns the Uniprot IDs which have missing diseases as well as the percentage of missing values, from the diseases dictionnary outputted by `get_uniprot_diseases`"""
    ids_missing_diseases = [
        k
        for k, v in diseases.items()
        if (v["comments"] == []) and (v["keywords"] == [])
    ]
    percentage_missing = len(ids_missing_diseases) / len(list(diseases.keys())) * 100
    return ids_missing_diseases, percentage_missing


def add_keywords_when_comments_missing(
    diseases: dict[str, dict[str, list[str]]]
) -> pd.DataFrame:
    df = pd.DataFrame(diseases).T.reset_index(
        names="UniProt (SwissProt) Primary ID of Target Chain"
    )
    df["comments"] = df["comments"].apply(lambda y: np.nan if len(y) == 0 else y)
    df["keywords"] = df["keywords"].apply(lambda y: np.nan if len(y) == 0 else y)
    df["comments_bfill"] = df[["comments", "keywords"]].bfill(axis=1)["comments"]
    return df


def weigh_each_comment(
    binding_db_df: pd.DataFrame, diseases_df: pd.DataFrame
) -> pd.DataFrame:
    """Weigh each comment by the number of occurences in the original binding db dataframe"""
    diseases_df = diseases_df.merge(
        binding_db_df["UniProt (SwissProt) Primary ID of Target Chain"]
        .value_counts()
        .reset_index()
    )
    diseases_df["comments_weighted"] = diseases_df["comments"] * diseases_df["count"]
    diseases_df["keywords_weighted"] = diseases_df["keywords"] * diseases_df["count"]
    diseases_df["comments_bfill_weighted"] = (
        diseases_df["comments_bfill"] * diseases_df["count"]
    )
    return diseases_df


def sort_diseases(disease_name: str) -> str:
    disease_name_lower = disease_name.lower()

    in_mapping = {
        "Neurodegeneration": [
            "parkinson",
            "alzeihmer",
        ],
        "Cancer": [
            "leukemia",
            "oncogene",
        ],
        "Obesity": [],
        "Neoplasia": [],
        "Hirschsprung Disease": [],
        "Immunodeficiency": [],
        "Agammaglobulinemia": [],
        "Thrombocythemia": [],
        "Diabetes": [],
        "Cornelia de Lange syndrome": [],
        "QT syndrome": [],
    }

    endswith_mapping = {"Cancer": ["oma"]}

    for k, v in in_mapping.items():
        if k.lower() in disease_name_lower:
            return k
        for v_i in v:
            if v_i in disease_name_lower:
                return k

    for k, v in endswith_mapping.items():
        for v_i in v:
            if disease_name_lower.endswith(v_i):
                return k

    return disease_name


def merge_and_explode_comments(
    binding_db_df: pd.DataFrame, diseases_df: pd.DataFrame
) -> pd.DataFrame:
    data = (
        binding_db_df.merge(diseases_df, how="inner")
        .dropna(subset="comments_bfill")
        .explode("comments_bfill")["comments_bfill"]
        .to_frame()
    )
    data["Disease Classes"] = data["comments_bfill"].apply(sort_diseases)
    return data
