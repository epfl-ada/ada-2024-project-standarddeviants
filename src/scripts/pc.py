import json

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from src.scripts import timeseries_anim
from src.scripts import smiles
from tqdm.notebook import tqdm


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

def prepare_pca(smiles_df, phase_4_count_df, citations_path, patents_path, uniprot_ids=["P07949","P14416"]):
    """Prepare data for PCA to overlay with success metrics, by filtering UniprotIDs and merging with citations and patents.

    Args:
        smiles_df (pd.DataFrame): BindingDB dataframe containing the following columns: "UniProt (SwissProt) Primary ID of Target Chain","Ligand SMILES","IC50 (nM)","Article DOI","Institution","Patent Number","ZINC ID of Ligand".
        phase_4_count_df (pd.DataFrame): contains the phase 4 and completed clinical trials per disease.
        citations_path (str): path to the citations.json file.
        patents_path (str): path to the patents.json file.
        uniprot_ids (list, str): List of UniprotIDs of interest. Defaults to ["P07949","P14416"].

    Returns:
        list: with two pd.DataFrame (result and targets).
    """
    targets = timeseries_anim.filter_uniprotID(smiles_df, uniprot_ids=uniprot_ids)
    targets = targets.drop_duplicates(subset=["Ligand SMILES"])
    targets = merge_citations_with_targets(targets, citations_path)
    targets = process_and_merge_patents(targets, patents_path)
    targets = pd.merge(
        left=targets,
        right=phase_4_count_df,
        left_on="ZINC ID of Ligand",
        right_on="ZINC ID of Ligand",
        how="outer",
    ).dropna(subset="Ligand SMILES")
    targets = clean_and_transform_IC50(targets, IC50_column="IC50 (nM)")
    targets = timeseries_anim.get_ligands_fingerprint(targets)
    store = []
    for target in targets["UniProt (SwissProt) Primary ID of Target Chain"].unique():
        temp = targets[
            targets["UniProt (SwissProt) Primary ID of Target Chain"] == target
        ].copy()
        temp, _, _ = timeseries_anim.PCA_fingerprints(temp)
        store.append(temp)
    result = pd.concat(store, ignore_index=True)
    result = result.drop_duplicates(subset=["Ligand SMILES"])
    return [result,targets]

def kmeans_selection(data, range):
    """Select the optimal number of clusters for KMeans by calculating silhouette scores and SSE for different k values.

    Args:
        data (pd.DataFrame): input data that will be used for clustering.
        range (iterable): integers representing the number of clusters to test.

    Returns:
        list: a list of dictionaries, each containing the k (number of clusters), corresponding silhouette score and sse.
    """
    scores = []
    for k in range:
        kmeans = KMeans(n_clusters=k, random_state=10).fit(data[["PC1", "PC2", "PC3"]])
        labels = kmeans.predict(data[["PC1", "PC2", "PC3"]])
        score = silhouette_score(data[["PC1", "PC2", "PC3"]], labels)
        scores.append({"k": k, "silhouette_score": score, "sse": kmeans.inertia_})
    return scores

def cluster_representatives(result_with_kmeans: pd.DataFrame):
    """Find the representative points for each cluster based on median value of the "log(IC50+1) (nM)" column.

    Args:
        result_with_kmeans (pd.DataFrame): contains cluster assignments and IC50.

    Returns:
        pd.DataFrame: contains the representative point for each cluster, sorted by log(IC50+1) values.
    """
    representatives = []
    for cluster in result_with_kmeans["cluster"].unique():
        cluster_results = result_with_kmeans.query(f"cluster == '{cluster}'").dropna(
            subset=["log(IC50+1) (nM)"]
        )
        cluster_results = cluster_results.sort_values("log(IC50+1) (nM)")
        representatives.append(cluster_results.iloc[len(cluster_results) // 2])
    return pd.DataFrame(representatives).sort_values("log(IC50+1) (nM)")

def get_representative_descriptors(representatives: pd.DataFrame):
    """Calculate molecular descriptors for each compound in the 'Ligand SMILES' column.

    Args:
        representatives (pd.DataFrame): input dataframe containing the 'Ligand SMILES' column.

    Returns:
        pd.DataFrame: contains the original data along with the calculated molecular descriptors.
    """
    tqdm.pandas()

    descriptors = representatives["Ligand SMILES"].progress_apply(
        smiles.get_MolDescriptors
    )
    isdict = lambda x: isinstance(x, dict)
    descriptors_df = pd.DataFrame(
        descriptors[descriptors.apply(isdict)].to_list(),
        index=descriptors[descriptors.apply(isdict)].index,
    )
    representatives_with_descriptors = pd.concat(
        [representatives, descriptors_df], axis=1
    )
    return representatives_with_descriptors