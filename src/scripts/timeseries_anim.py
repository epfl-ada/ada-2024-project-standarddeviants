import os

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from src.scripts import smiles
from src.utils.utils import group_categories


def filter_RET(
    smiles_df: pd.DataFrame, uniprot_ids: list[str] = ["P07949"]
) -> pd.DataFrame:
    return smiles_df[
        smiles_df["UniProt (SwissProt) Primary ID of Target Chain"].isin(uniprot_ids)
    ]


def get_ligands_fingerprint(df: pd.DataFrame) -> pd.DataFrame:
    target = df.copy(deep=True)
    target = target.dropna(subset=["Ligand SMILES"])
    target["Ligand Fingerprint"] = target["Ligand SMILES"].apply(smiles.get_fingerprint)
    target = target.dropna(subset=["Ligand Fingerprint"])
    return target


def merge_DOI_data(df: pd.DataFrame, doi_data_path: str = None) -> pd.DataFrame:
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
    target = df.copy(deep=True)
    data = np.stack(target["Ligand Fingerprint"].dropna())

    scaler = StandardScaler()
    data = scaler.fit_transform(data)

    pca = PCA(n_components=n_components, random_state=random_state)

    pca_coords = pca.fit_transform(data)
    target[["PC1", "PC2", "PC3"]] = pca_coords

    return target, scaler, pca


def group_institutions(df: pd.DataFrame) -> pd.DataFrame:
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
