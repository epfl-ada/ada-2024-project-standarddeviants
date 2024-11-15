import json
import os

import numpy as np
import pandas as pd

from src.utils.utils import group_categories


def load_uniprotid_diseases(
    filepath: str | os.PathLike | None = None,
) -> dict[str, dict[str, list[str]]]:
    """
    Load a dictionary of UniProt IDs and their associated diseases.
    
    Parameters:
        filepath (str | os.PathLike | None): Path to the JSON file containing the diseases dictionary.
                                              If None, defaults to 'data/UniprotID_disases.json'.
    
    Returns:
        dict[str, dict[str, list[str]]]: A dictionary where keys are UniProt IDs, and values contain 
                                         'comments' and 'keywords' lists describing diseases.
    """
    if filepath is None:
        filepath = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "data",
            "UniprotID_disases.json",
        )
    with open(filepath, "r") as f:
        diseases = json.load(f)
    return diseases


def quantify_missing_diseases(
    diseases: dict[str, dict[str, list[str]]]
) -> tuple[list[str], float]:
     """
    Identify UniProt IDs with missing disease data and calculate the percentage of missing values.
    
    Parameters:
        diseases (dict[str, dict[str, list[str]]]): The diseases dictionary with 'comments' and 'keywords'.
        
    Returns:
        tuple[list[str], float]: List of UniProt IDs missing diseases and the percentage of missing values.
    """
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
    """
    Create a DataFrame from the diseases dictionary, replacing missing comments with keywords if available.
    
    Parameters:
        diseases (dict[str, dict[str, list[str]]]): The diseases dictionary with 'comments' and 'keywords'.
        
    Returns:
        pd.DataFrame: A DataFrame with columns for 'comments', 'keywords', and a backfilled 'comments_bfill'.
    """
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
     """
    Weight each disease comment by its frequency in the binding database.
    
    Parameters:
        binding_db_df (pd.DataFrame): The binding database containing UniProt IDs.
        diseases_df (pd.DataFrame): DataFrame containing diseases with UniProt IDs, comments, and keywords.
        
    Returns:
        pd.DataFrame: A modified diseases DataFrame with weighted columns for 'comments', 'keywords', and 'comments_bfill'.
    """
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
    """
    Map a disease name to a general disease category based on specific keywords or patterns.
    
    Parameters:
        disease_name (str): The name of the disease to categorize.
        
    Returns:
        str: The category name for the disease.
    """
    in_mapping = {
        "Neurodegeneration": [
            "parkinson",
            "alzheimer",
        ],
        "Cancer": ["leukemia", "oncogene", "hemangioma", "glioma"],
        "Obesity": [],
        "Neoplasia": [],
        "Hirschsprung Disease": [],
        "Immunodeficiency": [],
        "Agammaglobulinemia": [],
        "Thrombocythemia": [],
        "Diabetes": [],
        "Cornelia de Lange syndrome": [],
        "QT syndrome": [],
        "Inflammatory skin and bowel disease": [],
        "Osteopetrosis": [],
        "Pregnancy loss": [],
        "Thrombophilia": [],
        "Periodic fever": [],
        "Lipodystrophy": [],
        "Schizophrenia": [],
    }

    endswith_mapping = {"Cancer": ["oma"]}

    return group_categories(
        disease_name, in_mapping=in_mapping, endswith_mapping=endswith_mapping
    )


def merge_and_explode_comments(
    
    binding_db_df: pd.DataFrame, diseases_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge the binding database with diseases, expand comments to individual rows, and categorize diseases.
    
    Parameters:
        binding_db_df (pd.DataFrame): The binding database DataFrame with UniProt IDs.
        diseases_df (pd.DataFrame): Diseases DataFrame with UniProt IDs and comments.
        
    Returns:
        pd.DataFrame: A DataFrame with exploded 'comments_bfill' column and a categorized 'Disease Classes' column.
    """
    data = (
        binding_db_df.merge(diseases_df, how="inner")
        .dropna(subset="comments_bfill")
        .explode("comments_bfill")["comments_bfill"]
        .to_frame()
    )
    data["Disease Classes"] = data["comments_bfill"].apply(sort_diseases)
    return data
