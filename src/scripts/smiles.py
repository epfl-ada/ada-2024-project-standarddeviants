# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import py3Dmol

# plotting tools
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Draw, rdCoordGen
from rdkit.Chem.Draw import IPythonConsole


def molecule_to_3d(molecule):
    mol = Chem.Mol(molecule)
    mol = AllChem.AddHs(mol, addCoords=True)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def get_Ki(df):
    """
    Extracts and converts Ki values from nanomolar (nM) format, removing any '>' or '<' symbols.

    Parameters:
        df (pd.DataFrame): DataFrame containing a column 'Ki (nM)' with inhibition constant values.

    Returns:
        pd.Series: Ki values as floats.
    """
    return (
        df["Ki (nM)"]
        .astype(str)
        .str.replace(">", "")
        .str.replace("<", "")
        .astype(float)
    )


def get_IC50(df):
    """
    Extracts and converts IC50 values from nanomolar (nM) format, removing any '>' or '<' symbols.

    Parameters:
        df (pd.DataFrame): DataFrame containing a column 'IC50 (nM)' with inhibition constant values.

    Returns:
        pd.Series: Ki values as floats.
    """
    return (
        df["IC50 (nM)"]
        .astype(str)
        .str.replace(">", "")
        .str.replace("<", "")
        .astype(float)
    )


##############  EXTRACTING MOLECULAR FEATURES FROM SMILES  ###############
# Fingerprint
def get_fingerprint(smiles: str):  # , radius:int=3, list_fmt=True)->list:
    """
    Generates the Morgan fingerprint for a molecule from its SMILES representation.

    Parameters:
        smiles (str): SMILES string of the molecule.

    Returns:
        rdkit.DataStructs.cDataStructs.ExplicitBitVect or NaN: Fingerprint object or NaN if molecule is invalid.
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
    """
    Calculates the number of hydrogen bond donors in a molecule from its SMILES.

    Parameters:
        smiles (str): SMILES string of the molecule.

    Returns:
        int or False: Number of hydrogen bond donors or False if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return Descriptors.NumHDonors(mol)


def get_Hacceptors(smiles):
    """
    Calculates the number of hydrogen bond acceptors in a molecule from its SMILES.

    Parameters:
        smiles (str): SMILES string of the molecule.

    Returns:
        int or False: Number of hydrogen bond acceptors or False if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return Descriptors.NumHAcceptors(mol)


def get_MW(smiles):
    """
    Calculates the molecular weight of a molecule from its SMILES.

    Parameters:
        smiles (str): SMILES string of the molecule.

    Returns:
        float or False: Molecular weight or False if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return Descriptors.MolWt(mol)


def get_LogP(smiles):
    """
    Calculates the LogP (partition coefficient) of a molecule from its SMILES.

    Parameters:
        smiles (str): SMILES string of the molecule.

    Returns:
        float or False: LogP value or False if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return Descriptors.MolLogP(mol)


def get_MolDescriptors(smiles: str):
    """
    Calculates the LogP (partition coefficient) of a molecule from its SMILES.

    Parameters:
        smiles (str): SMILES string of the molecule.

    Returns:
        dict or False: Dictionary of all molecular descriptors or False if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "Molecular Weight": None,
            "C LogP": None,
            "FractionCSP3": None,
            "HallKierAlpha": None,
            "HeavyAtomCount": None,
            "HeavyAtomMolWt": None,
            "Ipc": None,
            "Kappa1": None,
            "Kappa2": None,
            "Kappa3": None,
            "LabuteASA": None,
            "MaxAbsEStateIndex": None,
            "MaxAbsPartialCharge": None,
            "MaxEStateIndex": None,
            "MaxPartialCharge": None,
            "MinAbsEStateIndex": None,
            "MinAbsPartialCharge": None,
            "MinEStateIndex": None,
            "MinPartialCharge": None,
            "MolLogP": None,
            "MolMR": None,
            "NHOHCount": None,
            "NOCount": None,
            "NumAliphaticCarbocycles": None,
            "NumAliphaticHeterocycles": None,
            "NumAliphaticRings": None,
            "NumAromaticCarbocycles": None,
            "NumAromaticHeterocycles": None,
            "NumAromaticRings": None,
            "NumHAcceptors": None,
            "NumHDonors": None,
            "NumHeteroatoms": None,
            "NumRadicalElectrons": None,
            "NumRotatableBonds": None,
            "NumSaturatedCarbocycles": None,
            "NumSaturatedHeterocycles": None,
            "NumSaturatedRings": None,
            "NumValenceElectrons": None,
        }
    else:
        return {
            "Molecular Weight": Descriptors.MolWt(mol),
            "C LogP": Descriptors.MolLogP(mol),
            "FractionCSP3": Descriptors.FractionCSP3(mol),
            "HallKierAlpha": Descriptors.HallKierAlpha(mol),
            "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),
            "HeavyAtomMolWt": Descriptors.HeavyAtomMolWt(mol),
            "Ipc": Descriptors.Ipc(mol),
            "Kappa1": Descriptors.Kappa1(mol),
            "Kappa2": Descriptors.Kappa2(mol),
            "Kappa3": Descriptors.Kappa3(mol),
            "LabuteASA": Descriptors.LabuteASA(mol),
            "MaxAbsEStateIndex": Descriptors.MaxAbsEStateIndex(mol),
            "MaxAbsPartialCharge": Descriptors.MaxAbsPartialCharge(mol),
            "MaxEStateIndex": Descriptors.MaxEStateIndex(mol),
            "MaxPartialCharge": Descriptors.MaxPartialCharge(mol),
            "MinAbsEStateIndex": Descriptors.MinAbsEStateIndex(mol),
            "MinAbsPartialCharge": Descriptors.MinAbsPartialCharge(mol),
            "MinEStateIndex": Descriptors.MinEStateIndex(mol),
            "MinPartialCharge": Descriptors.MinPartialCharge(mol),
            "MolLogP": Descriptors.MolLogP(mol),
            "MolMR": Descriptors.MolMR(mol),
            "NHOHCount": Descriptors.NHOHCount(mol),
            "NOCount": Descriptors.NOCount(mol),
            "NumAliphaticCarbocycles": Descriptors.NumAliphaticCarbocycles(mol),
            "NumAliphaticHeterocycles": Descriptors.NumAliphaticHeterocycles(mol),
            "NumAliphaticRings": Descriptors.NumAliphaticRings(mol),
            "NumAromaticCarbocycles": Descriptors.NumAromaticCarbocycles(mol),
            "NumAromaticHeterocycles": Descriptors.NumAromaticHeterocycles(mol),
            "NumAromaticRings": Descriptors.NumAromaticRings(mol),
            "NumHAcceptors": Descriptors.NumHAcceptors(mol),
            "NumHDonors": Descriptors.NumHDonors(mol),
            "NumHeteroatoms": Descriptors.NumHeteroatoms(mol),
            "NumRadicalElectrons": Descriptors.NumRadicalElectrons(mol),
            "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
            "NumSaturatedCarbocycles": Descriptors.NumSaturatedCarbocycles(mol),
            "NumSaturatedHeterocycles": Descriptors.NumSaturatedHeterocycles(mol),
            "NumSaturatedRings": Descriptors.NumSaturatedRings(mol),
            "NumValenceElectrons": Descriptors.NumValenceElectrons(mol),
        }


def tanimoto(fp1, fp2) -> float:
    """
    Computes the Tanimoto similarity between two molecular fingerprints.

    Parameters:
        fp1, fp2: Molecular fingerprints.

    Returns:
        float: Tanimoto similarity score between the fingerprints.
    """
    return DataStructs.FingerprintSimilarity(fp1, fp2)


def tanimoto_smiles(sm1, sm2) -> float:
    """
    Computes the Tanimoto similarity between two molecules based on their SMILES strings.

    Parameters:
        sm1, sm2 (str): SMILES strings of the molecules.

    Returns:
        float: Tanimoto similarity score.
    """
    return DataStructs.FingerprintSimilarity(get_fingerprint(sm1), get_fingerprint(sm2))


def tanimoto_matrix(df: pd.DataFrame, index: str = "Ligand SMILES") -> pd.DataFrame:
    """
    Creates a Tanimoto similarity matrix for a DataFrame of molecules based on their SMILES.

    Parameters:
        df (pd.DataFrame): DataFrame with a column of SMILES strings.
        index (str): Column name for SMILES strings (default is 'Ligand SMILES').

    Returns:
        pd.DataFrame: Square DataFrame of Tanimoto similarity scores.
    """
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
    """
    Finds two closest divisors of n to create a grid layout.

    Parameters:
        n (int): Number to find divisors for.

    Returns:
        tuple: Closest divisors for n as (rows, columns).
    """
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
    """
    Plots a heatmap of the Tanimoto similarity matrix.

    Parameters:
        tanimoto_matrix (pd.DataFrame): Tanimoto similarity matrix.
        show_every_n_labels (int): Interval for showing labels.
        label_type (str): Label type for the plot.
        title (str): Title of the plot.

    Returns:
        pd.DataFrame: Displayed Tanimoto similarity matrix.
    """
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
    """
    Displays molecular structures with corresponding target names.

    Parameters:
        df (pd.DataFrame): DataFrame with 'Ligand SMILES' and 'Target Name' columns.
        title (str): Title for the display.
        nrow (int): Number of rows in the grid.
        ncol (int): Number of columns in the grid.
        grouped_by (int): Number of images grouped by target.
    """
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
    """
    Plots the most potent ligands for each target.

    Parameters:
        df_with_unique_targets (pd.DataFrame): DataFrame with unique target names.
        df_with_chems (pd.DataFrame): DataFrame with chemical data.
        chems_per_target (int): Number of ligands per target to plot.
        title (str): Title for the plot.
        nrow (int): Number of rows in the plot grid.
        ncol (int): Number of columns in the plot grid.
    """
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
    """
    Displays molecular structures in a grid.

    Parameters:
        df (pd.DataFrame): DataFrame with 'Ligand SMILES' column.
        title (str): Title for the plot.
        n_rows (int): Number of rows in the grid.
        n_cols (int): Number of columns in the grid.
        random_sample (bool): Whether to randomly sample SMILES strings.
    """

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


def lipinski(smiles: str, verbose=False) -> bool:
    """
    Checks if a molecule adheres to Lipinski's Rule of Five.

    Parameters:
        smiles (str): SMILES string of the molecule.
        verbose (bool): Whether to print details of the rules.

    Returns:
        bool: True if molecule follows Lipinski's rules, False otherwise.
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


# 3D stuff
def draw3D(smiles, outpath=None, display=False):
    mol = Chem.MolFromSmiles(smiles)
    rdCoordGen.AddCoords(mol)
    mol = AllChem.AddHs(mol, addCoords=True)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    view = py3Dmol.view(
        data=Chem.MolToMolBlock(mol),  # Convert the RDKit molecule for py3Dmol
        style={"stick": {}, "sphere": {"scale": 0.3}},
    )
    view.setBackgroundColor("#222529")
    if display:
        view.zoomTo()
        view.show()
    if outpath is not None:
        view.write_html(f=outpath, fullpage=True)
    return view.write_html(f=outpath, fullpage=True)
