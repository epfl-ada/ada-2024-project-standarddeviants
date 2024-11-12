# imports
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# plotting tools
import seaborn as sns
import umap
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Draw

# rdkit tools
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA

# import networkx as nx


def get_Ki(df):
    return (
        df["Ki (nM)"]
        .astype(str)
        .str.replace(">", "")
        .str.replace("<", "")
        .astype(float)
    )


##############  EXTRACTING MOLECULAR FEATURES FROM SMILES  ###############
# Fingerprint
def get_fingerprint(smiles: str):  # , radius:int=3, list_fmt=True)->list:
    """Gets the fingerprint from a SMILES string. Easy to use with pd.Series.apply() on a column with smiles data.

    Args:
        smiles (str): SMILES string of molecule of interest
        radius (int, optional): Radius of the generated fingerprint. Defaults to 3.

    Returns:
        list: Fingerprint of molecule of interest, Nan of no molecule or fingerprint is found.
    """
    list_fmt = False
    radius = 3
    fpgen = AllChem.GetMorganGenerator(radius=radius)
    mol = Chem.MolFromSmiles(smiles)
    if not mol is None:
        if list_fmt:
            return fpgen.GetFingerprint(mol).ToList()
        else:
            return fpgen.GetFingerprint(mol)
    else:
        return np.nan


def get_Hdonors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return Descriptors.NumHDonors(mol)


def get_Hacceptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return Descriptors.NumHAcceptors(mol)


def get_MW(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return Descriptors.MolWt(mol)


def get_LogP(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return Descriptors.MolLogP(mol)


# Compute the Tanimoto similarity between all pairs
def tanimoto(fp1, fp2) -> float:
    return DataStructs.FingerprintSimilarity(fp1, fp2)


def tanimoto_smiles(sm1, sm2) -> float:
    return DataStructs.FingerprintSimilarity(get_fingerprint(sm1), get_fingerprint(sm2))


def tanimoto_matrix(df: pd.DataFrame, index: str = "Ligand SMILES") -> pd.DataFrame:
    # fingerprints
    fingerprints = df["Ligand SMILES"].apply(get_fingerprint).dropna()
    # Compute the Tanimoto similarity between all pairs
    return pd.DataFrame(
        np.array(
            [[tanimoto(fp1, fp2) for fp2 in fingerprints] for fp1 in fingerprints]
        ),
        index=df[index],
        columns=df[index],
    )


# PLOTTING FUNCTIONS
def closest_dividers(n: int) -> tuple:
    import math

    dividers = []
    for i in range(1, n + 1):
        if n % i == 0:
            dividers.append(i)
    mid = int(len(dividers) / 2)
    if math.sqrt(n).is_integer():
        return (dividers[mid], dividers[mid])
    else:
        return (dividers[mid - 1], dividers[mid])


def tanimoto_plot(
    tanimoto_matrix,
    show_every_n_labels: int = 0,
    label_type: str = "Target Name",
    title: str = "Tanimoto Similarity Matrix",
):
    fig, ax = plt.subplots(1, 1)

    # plot labels
    if show_every_n_labels:
        new_labels = []
        for idx, name in enumerate(chem_df[label_type]):
            if idx % show_every_n_labels == 0:
                new_labels.append(name)
            else:
                new_labels.append("")
    else:
        new_labels = ""

    # plot
    sns.heatmap(
        tanimoto_matrix,
        cmap="bone",
        xticklabels=new_labels,
        yticklabels=new_labels,
        cbar_kws={"location": "left"},
        ax=ax,
    )
    ax.set_title(title)
    ax.yaxis.tick_right()
    plt.yticks(rotation=0)
    plt.show()

    return tanimoto_matrix


def show_smiles_with_target(
    df: pd.DataFrame,
    title: str = "",
    nrow: int = None,
    ncol: int = None,
    grouped_by: int = None,
):
    smiles_and_target = df[["Ligand SMILES", "Target Name"]].reset_index()
    length = smiles_and_target.shape[0]
    if nrow and ncol:
        n_rows = nrow
        n_cols = ncol
    else:
        n_cols, n_rows = closest_dividers(length)

    fig, axes = plt.subplots(
        nrows=n_rows, ncols=n_cols, figsize=(2 * n_cols, 2 * n_rows)
    )

    for ax, i in zip(axes.flat, np.arange(length)):
        mol = Chem.MolFromSmiles(smiles_and_target.loc[i, "Ligand SMILES"])
        tar = smiles_and_target.loc[i, "Target Name"]
        img = Draw.MolToImage(mol)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(f"Target: {tar}")

    if grouped_by:
        for ax, i in zip(axes.flat, np.arange(length)):
            if (i - 1) % grouped_by != 0:
                ax.set_title("")

    fig.suptitle(title, fontweight="bold")
    plt.tight_layout()
    plt.show()


def plot_potent_chems(
    df_with_unique_targets: pd.DataFrame,
    df_with_chems: pd.DataFrame,
    chems_per_target: int = 1,
    title: str = "Most Potent Ligand(s) of Each Target",
    nrow: int = None,
    ncol: int = None,
) -> None:
    most_potent_stored = []  # store all structures
    for target in df_with_unique_targets.index:
        most_potent = df_with_chems[df_with_chems["Target Name"] == target].sort_values(
            by="IC50 (nM)", ascending=True
        )[:chems_per_target]
        most_potent_stored.append(most_potent)

    most_potent_stored = pd.concat(most_potent_stored)[
        ["Ligand SMILES", "Target Name", "IC50 (nM)"]
    ]
    show_smiles_with_target(
        most_potent_stored,
        title=title,
        grouped_by=chems_per_target,
        ncol=ncol,
        nrow=nrow,
    )


def show_smiles(df: pd.DataFrame, title="", n_rows=5, n_cols=5, random_sample=False):
    unique_smiles = df["Ligand SMILES"].unique()

    n_sampled = n_rows * n_cols
    if random_sample:
        smiles = np.random.choice(unique_smiles, n_sampled, replace=False)
    else:
        smiles = unique_smiles[:n_sampled]

    fig, axes = plt.subplots(
        nrows=n_rows, ncols=n_cols, figsize=(2 * n_cols, 2 * n_rows)
    )

    for ax, i in zip(axes.flat, np.arange(n_sampled)):
        mol = Chem.MolFromSmiles(smiles[i])
        img = Draw.MolToImage(mol)
        ax.imshow(img)
        ax.axis("off")
    fig.suptitle(title)
    plt.tight_layout()
    plt.show()


def respects_lipinski(smiles: str, verbose=False) -> bool:
    """
    Checks if a molecule respects Lipinski's Rule of Five.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule respects Lipinski's rules, False otherwise.
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES string.")
        return False

    # Calculate Lipinski properties
    hbd = Descriptors.NumHDonors(mol)  # Hydrogen bond donors
    hba = Descriptors.NumHAcceptors(mol)  # Hydrogen bond acceptors
    mw = Descriptors.MolWt(mol)  # Molecular weight
    logp = Descriptors.MolLogP(mol)  # Partition coefficient (LogP)

    # Check Lipinski’s rules
    respects_rules = (
        hbd <= 5
        and hba <= 10  # No more than 5 H-bond donors
        and mw < 500  # No more than 10 H-bond acceptors
        and logp <= 5  # Molecular weight less than 500 Da  # LogP not greater than 5
    )

    # Print details
    if verbose:
        print(f"SMILES: {smiles}")
        print(f"Hydrogen Bond Donors (HBD): {hbd} (<= 5)")
        print(f"Hydrogen Bond Acceptors (HBA): {hba} (<= 10)")
        print(f"Molecular Weight (MW): {mw:.2f} (< 500 Da)")
        print(f"LogP: {logp:.2f} (<= 5)")
        print(f"Respects Lipinski's Rule of Five: {respects_rules}")
    return respects_rules
