import json

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def merge_citations_with_targets(targets, citations_path="../src/data/citations.json"):
    """
    Merges the targets DataFrame with a citations JSON file.

    Args:
        targets (pd.DataFrame): The targets DataFrame to merge with citations.
        citations_path (str): Path to the citations JSON file.

    Returns:
        pd.DataFrame: The merged and transformed DataFrame.
    """
    # Load the citations JSON
    with open(citations_path, "r") as f:
        citations = json.load(f)
    citations = pd.DataFrame(citations)

    # Merge and transform
    merged_df = (
        pd.merge(
            left=targets,
            right=citations,
            left_on="Article DOI",
            right_on="doi",
            how="outer",
        )
        .dropna(subset=["Ligand SMILES"])
        .drop(columns=["Article DOI"])
    )

    return merged_df


def process_and_merge_patents(targets, patents_path="../src/data/patents.json"):
    """
    Processes a patents JSON file into a DataFrame and merges it with the targets DataFrame.
    Args:
        targets (pd.DataFrame): The targets DataFrame to merge with the patents data.
        patents_path (str): Path to the patents JSON file.

    Returns:
        pd.DataFrame: The merged and transformed DataFrame.
    """
    # patents
    with open(patents_path, "r") as f:
        patents = json.load(f)
    patents_df = pd.DataFrame(
        [
            {
                "patent": patent["patent"],
                "patent_status": patent["info"].get("status", np.nan),
                "families citing": int(patent["info"].get("families citing", 0) or 0),
                "cited by": int(patent["info"].get("cited by", 0) or 0),
            }
            for patent in patents
            if isinstance(patent, dict) and isinstance(patent.get("info"), dict)
        ]
    )

    patents_df["patent_citations"] = (
        patents_df["families citing"] + patents_df["cited by"]
    )
    patents_df.drop(columns=["families citing", "cited by"], inplace=True)

    # Merge
    merged_df = (
        pd.merge(
            left=targets,
            right=patents_df,
            left_on="Patent Number",
            right_on="patent",
            how="outer",
        )
        .dropna(subset=["Ligand SMILES"])
        .drop(columns=["Patent Number"])
    )

    return merged_df


def clean_and_transform_IC50(
    df: pd.DataFrame,
    IC50_column: str = "IC50 (nM)",
    log_column_name: str = "log(IC50+1) (nM)",
) -> pd.DataFrame:
    """
    Cleans the IC50 column in the DataFrame and calculates the log(IC50+1).

    Args:
        df (pd.DataFrame): The input DataFrame containing IC50 data.
        IC50_column (str): The name of the IC50 column to clean.
        log_column_name (str): The name of the column to store the log-transformed IC50.

    Returns:
        pd.DataFrame: The DataFrame with cleaned IC50 values and the log-transformed column.
    """
    # Handle string nan values
    df.replace(" NV,", np.nan, inplace=True)

    # Clean IC50 data
    df[IC50_column] = df[IC50_column].astype(str).str.replace(" C", "")
    df[IC50_column] = (
        df[IC50_column].str.replace(">", "").str.replace("<", "").astype(float)
    )

    # Add column
    df[log_column_name] = (df[IC50_column] + 1).apply(np.log10)

    return df


def kmeans_selection(data, range):
    """
    Evaluates KMeans clustering for a range of cluster numbers using silhouette scores and sum of squared errors (SSE).

    Args:
        data (pd.DataFrame): DataFrame with "PC1", "PC2", and "PC3" columns for clustering.
        range (iterable): Range of cluster counts (k) to evaluate.

    Returns:
        list: A list of dictionaries with:
            - "k": Number of clusters.
            - "silhouette_score": Silhouette score for k clusters.
            - "sse": Sum of squared errors (SSE) for k clusters.
    """
    scores = []
    for k in range:
        kmeans = KMeans(n_clusters=k, random_state=10).fit(data[["PC1", "PC2", "PC3"]])
        labels = kmeans.predict(data[["PC1", "PC2", "PC3"]])
        score = silhouette_score(data[["PC1", "PC2", "PC3"]], labels)
        scores.append({"k": k, "silhouette_score": score, "sse": kmeans.inertia_})
    return scores
