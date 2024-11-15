from ete3 import NCBITaxa

from src.utils import utils

# Taxonomy
ncbi = NCBITaxa()


def get_taxon_id(species_name):
    """
    Get the NCBI taxon ID for a given species name using ete3.

    Parameters:
        species_name (str): The scientific name of the species.

    Returns:
        taxon_id (int): The corresponding taxon ID, or None if not found.
    """
    try:
        # Get taxon ID for the species name
        taxon_id = ncbi.get_name_translator([species_name])

        if species_name in taxon_id:
            return taxon_id[species_name][0]  # Return the taxon ID
        else:
            return None  # Species not found
    except Exception as e:
        print(f"Error fetching taxon ID: {e}")
        return None


# Tools for sequence data
def get_names(df):
    """
    Retrieve the recommended or submitted name of target chains.

    Parameters:
        df (pd.DataFrame): DataFrame containing names of proteins in columns ['UniProt (SwissProt) Recommended Name of Target Chain',
                            'UniProt (TrEMBL) Submitted Name of Target Chain']
    Returns:
        pd.Series: Target chain names.
    """
    return df["UniProt (SwissProt) Recommended Name of Target Chain"].fillna(
        df["UniProt (TrEMBL) Submitted Name of Target Chain"]
    )


def get_id(df):
    """
    Retrieve the primary ID of target chains.

    Parameters:
        df (pd.DataFrame): DataFrame containing sequence data.

    Returns:
        pd.Series: Target chain primary IDs.
    """
    return df["UniProt (SwissProt) Primary ID of Target Chain"].fillna(
        df["UniProt (TrEMBL) Primary ID of Target Chain"]
    )


def to_fasta(df, filepath):
    """
    Save sequences as a FASTA file.

    Parameters:
        df (pd.DataFrame): DataFrame with sequence data.
        filepath (str): Path to save the FASTA file.
    """
    with open(filepath, "w") as file:
        for name, id, seq in zip(
            get_names(df), get_id(df), df["BindingDB Target Chain Sequence"]
        ):
            file.write(f">{id}|{name}\n{seq}\n")


def get_target_class(names_df):
    """
    Retrieve the class of a target protein.

    Parameters:
        df (pd.DataFrame): DataFrame containing names of proteins in columns ['UniProt (SwissProt) Recommended Name of Target Chain',
                            'UniProt (TrEMBL) Submitted Name of Target Chain']

    Returns:
        pd.Series: Target classes after grouping
    """
    mapping = {
        "Growth Factor Receptor": ["gfr"],
        "CDK": ["cyclin-dependent kinase"],
        "Polyprotein": [],
        "Histone Modifier": [
            "histone deacetylase",
            "bromodomain-containing protein",
            "histone demethylase",
        ],
        "RTK": [
            "tyrosine-protein kinase receptor",
            "tyrosine kinase receptor",
            "receptor tyrosine kinase",
        ],
        "Neurotransmitter receptor": [
            "serotonin receptor",
            "dopamine receptor",
            "hydroxytryptamine receptor",
            "dopamine receptor",
            "cannabinoid receptor",
            "glutamate receptor",
        ],
        "Neurotransmitter transporter": [
            "serotonin transporter",
            "dopamine transporter",
            "acetylcholinesterase",
            "cholinesterase",
        ],
        "Neurotransmitter synthesis": ["cholinesterase", "amine oxidase"],
        "Neural Peptide Receptor": ["orexin", "hypocretin"],
        "Hormone Receptor": [
            "corticotropin-releasing factor receptor",
            "prostaglandin",
            "estrogen receptor",
            "glucocorticoid receptor",
            "insulin receptor",
            "glucagon receptor",
            "vasopressin v1a receptor",
        ],
        "Caspase": [],
        "Blood Homeostasis": ["prothrombin", "coagulation factor", "plasma kallikrein"],
        "Non Receptor Tyr Kinase": ["btk", "tyk2", "jak", "syk", "abl1", "b-raf"],
        "Other Protein Kinase": [
            "glycogen synthase kinase-3",
            "rho-associated protein kinase 2",
            "raf proto-oncogene serine/threonine-protein kinase",
            "leucine-rich repeat serine/threonine-protein kinase 2",
            "serine/threonine-protein kinase pim-1",
            "receptor-interacting serine/threonine-protein kinase 1",
            "mitogen-activated protein kinase",
        ],
        "Transcription Factor": [
            "nuclear receptor ror-gamma",
            "peroxisome proliferator-activated receptor gamma",
        ],
        "Opioid Receptor": [],
        "Ion Channel": [
            "voltage-gated channel",
            "sodium channel",
            "potassium channel",
            "cation channel",
            "anion channel",
            "calcium channel",
        ],
        "Toll-Like Receptor": [],
        "Purine Receptor": ["adenosine receptor", "purinoceptor"],
        "Beta-Secretase": [],
        "Carbonic Anhydrase": [],
        "Protease": [],
        "_Cyclin_": ["cyclin-k", "cyclin-d1", "cyclin-t1", "cyclin-h", "cyclin-c"],
        "Ubiquitin Ligase": ["ubiquitin-protein ligase"],
        "Phosphatidylinositol Kinase": ["phosphatidylinositol"],
        "Interleukin Receptor": ["interleukin-1 receptor", "interleukin"],
        "NADP Regeneration": ["isocitrate dehydrogenase"],
        "Apoptosis Regulation": [
            "induced myeloid leukemia cell differentiation protein",
            "mcl-1",
        ],
        "GTPase": [],
        "Phosphodiesterase": [],
        "Integrase": [],
        "Complement Factor": [],
        "H3 Receptor": [],
        "Oxidase": ["cytochrome p450 3a4"],
    }
    mapper = lambda x: utils.group_categories(x, in_mapping=mapping)
    return get_names(names_df).dropna().astype(str).apply(mapper)
