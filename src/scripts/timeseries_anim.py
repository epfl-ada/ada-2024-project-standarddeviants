import os

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from src.scripts import smiles
from src.utils.utils import group_categories


def filter_uniprotID(
    smiles_df: pd.DataFrame, uniprot_ids: list[str] = ["P07949"]
) -> pd.DataFrame:
    """
    Filters a DataFrame to include rows with specified UniProt IDs.

    Args:
        smiles_df (pd.DataFrame): The input DataFrame containing UniProt IDs.
        uniprot_ids (list[str]): List of UniProt IDs to filter by (default: ["P07949"]).

    Returns:
        pd.DataFrame: Filtered DataFrame with rows matching the specified UniProt IDs.
    """
    return smiles_df[
        smiles_df["UniProt (SwissProt) Primary ID of Target Chain"].isin(uniprot_ids)
    ]


def get_ligands_fingerprint(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generates ligand fingerprints for rows with valid SMILES strings.

    Args:
        df (pd.DataFrame): Input DataFrame containing a "Ligand SMILES" column.

    Returns:
        target (pd.DataFrame): A DataFrame with a new "Ligand Fingerprint" column, filtered to include only rows with valid fingerprints.
    """
    target = df.copy(deep=True)
    target = target.dropna(subset=["Ligand SMILES"])
    target["Ligand Fingerprint"] = target["Ligand SMILES"].apply(smiles.get_fingerprint)
    target = target.dropna(subset=["Ligand Fingerprint"])
    return target


def merge_DOI_data(df: pd.DataFrame, doi_data_path: str = None) -> pd.DataFrame:
    """
    Merges a DataFrame with DOI metadata from a CSV file.

    Args:
        df (pd.DataFrame): Input DataFrame containing an "Article DOI" column.
        doi_data_path (str): Path to the CSV file with DOI metadata.

    Returns:
        pd.DataFrame: A DataFrame with merged DOI metadata, filtered for rows with valid DOIs.
    """
    if doi_data_path is None:
        doi_data_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "data",
            "metadata.csv",
        )
    target = df.copy(deep=True).dropna(subset=["Article DOI"])
    doi_metadata = pd.read_csv(doi_data_path).dropna()
    target = target.merge(doi_metadata)
    return target


def filter_institutions(df: pd.DataFrame, num_diff_years: int = 2) -> pd.DataFrame:
    """
    Filters a DataFrame to include only roes from institution present in at least a specifief nunmber of years

    Args :
        df (pd.DataFrame): Input Daframe containing an "Institution" column.
        num_diff_years (int) : Minimum number of unique years an institution must appear in.

    Returns :
        target (pd.DataFrame): Filtered DataFrame with rows from valid institutions.
    """
    target = df.copy(deep=True)
    # Count unique years for each institution
    institution_year_counts = target.groupby("Institution")["year"].nunique()
    # Filter institutions with at least 2 different years
    valid_institutions = institution_year_counts[
        institution_year_counts >= num_diff_years
    ].index
    # Filter the dataframe
    target = target[target["Institution"].isin(valid_institutions)]
    return target


def PCA_fingerprints(
    df: pd.DataFrame, n_components: int = 3, random_state: int = 0
) -> tuple[pd.DataFrame, StandardScaler, PCA]:
    """
    Applies PCA to ligand fingerprints and adds the PCA components as columns to the DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame containing a "Ligand Fingerprint" column.
        n_components (int): Number of PCA components to compute (default: 3).
        random_state (int): Random state for reproducibility (default: 0).

    Returns:
        target (pd.DataFrame): DataFrame with PCA components added as "PC1", "PC2", "PC3".
        scaler (StandardScaler): The scaler used to normalize the fingerprint data.
        pca (PCA): The PCA object used to compute the principal components.
    """
    target = df.copy(deep=True)
    data = np.stack(target["Ligand Fingerprint"].dropna())

    scaler = StandardScaler()
    data = scaler.fit_transform(data)

    pca = PCA(n_components=n_components, random_state=random_state)

    pca_coords = pca.fit_transform(data)
    target[["PC1", "PC2", "PC3"]] = pca_coords

    return target, scaler, pca


def group_institutions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Groups institutions into predefined categories based on keyword mappings and
    cleans the "Institution" column.

    Args:
        df (pd.DataFrame): Input DataFrame containing an "Institution" column.

    Returns:
        target (pd.DataFrame): DataFrame with the "Institution" column updated to reflect grouped categories, unmatched entries filled with NaN.
    """

    target = df.copy(deep=True)
    in_mapping = {
        "Pfizer": [],
        "MSD": ["Dohme"],
        "Bristol-Myers Squibb": [],
        "Amgen": [],
        "Novartis": [],
        "Janssen": [],
        "Eli Lilly": ["lilly"],
        "Roche": [],
        "Incyte": [],
        "Gilead": [],
        "Bayer": [],
        "Abbott": [],
        "Scripps Research Institute": ["scripps"],
        "The Burnham Institute": ["burnham"],
        "Genentech": [],
        "GlaxoSmithKline": ["gsk"],
        "Astrazeneca": [],
        "Abbvie": [],
        "Merck": [],
        "Boehring": [],
    }

    f = lambda x: group_categories(
        str(x),
        in_mapping=in_mapping,
        check_key_for_in_mapping=True,
    )

    target["Institution"] = (
        target["Institution"].apply(f).replace("TBA", np.nan).replace("nan", np.nan)
    )  # tba = to be attributed
    return target


def get_cumulative_years_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generates a cumulative DataFrame where rows are repeated for each year up to and including the row's year.

    Args:
        df (pd.DataFrame): Input DataFrame containing a "year" column.

    Returns:
        pd.DataFrame: A cumulative DataFrame with a "current_year" column added, containing repeated rows for all years up to the current year.
    """
    target = df.copy(deep=True)
    target = target.sort_values(by="year")

    return pd.concat(
        [
            target[target["year"] <= year].assign(
                current_year=year
            )  # Assign a new column called current year
            for year in target["year"].unique()
        ]
    )


def include_dummy_row_for_anim(
    df: pd.DataFrame, color_column: str, animation_frame: str = "current_year"
) -> pd.DataFrame:
    """
    Fills missing PCA coordinates with Nan, points won't be plotted but colors will be allocated

    Args:
        df (pd.DataFrame): Input DataFrame containing PCA and animation data.
        color_column (str): The column representing the feature.
        animation_frame (str): The column representing the animation.
    Returns:
        df_ (pd.DataFrame): A DataFrame with dummy rows added for missing combinations of colors and animation frames, with NaN for PCA coordinates.
    """
    df_ = df.copy(deep=True)

    # Get all unique institutions
    unique_colors = df_[color_column].unique()

    # Create a complete DataFrame of all institutions and years
    frames = df_[animation_frame].unique()
    color_frame_grid = pd.DataFrame(
        [(color, frame) for color in unique_colors for frame in frames],
        columns=[color_column, animation_frame],
    )

    # Merge with target_cumulative to ensure all combinations exist
    df_ = color_frame_grid.merge(df_, on=[color_column, animation_frame], how="left")

    # Fill missing PCA coordinates with NaN (points won't be plotted but colors will be allocated)
    df_[["PC1", "PC2", "PC3"]] = df_[["PC1", "PC2", "PC3"]].fillna(value=np.nan)

    return df_


def clean_IC50(df: pd.DataFrame, IC50_column: str = "IC50 (nM)") -> pd.DataFrame:
    """
    Cleans IC50 data by handling missing values, removing unwanted characters, and
    converting the column to float.

    Args:
        df (pd.DataFrame): Input DataFrame containing IC50 data.
        IC50_column (str): The name of the IC50 column to clean (default: "IC50 (nM)").

    Returns:
        pd.DataFrame: A DataFrame with cleaned IC50 data.
    """

    target_cumulative_with_IC50 = df.dropna(subset=[IC50_column])
    # clean binding kinetics data
    target_cumulative_with_IC50.replace(" NV,", np.nan, inplace=True)
    for col in [IC50_column]:
        target_cumulative_with_IC50[col] = (
            target_cumulative_with_IC50[col].astype(str).str.replace(" C", "")
        )
        target_cumulative_with_IC50[col] = (
            target_cumulative_with_IC50[col]
            .astype(str)
            .str.replace(">", "")
            .str.replace("<", "")
            .astype(float)
        )
    return target_cumulative_with_IC50
