from ete3 import NCBITaxa



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
def get_name(df):
    """
    Retrieve the recommended or submitted name of target chains.

    Parameters:
        df (pd.DataFrame): DataFrame containing sequence data.

    Returns:
        pd.Series: Target chain names.
    """
    return df['UniProt (SwissProt) Recommended Name of Target Chain'].fillna(df['UniProt (TrEMBL) Submitted Name of Target Chain'])

def get_id(df):
    """
    Retrieve the primary ID of target chains.

    Parameters:
        df (pd.DataFrame): DataFrame containing sequence data.

    Returns:
        pd.Series: Target chain primary IDs.
    """
    return df['UniProt (SwissProt) Primary ID of Target Chain'].fillna(df['UniProt (TrEMBL) Primary ID of Target Chain'])

def to_fasta(df, filepath):
    """
    Save sequences as a FASTA file.

    Parameters:
        df (pd.DataFrame): DataFrame with sequence data.
        filepath (str): Path to save the FASTA file.
    """
    with open(filepath, 'w') as file:
        for name, id, seq in zip(get_name(df), get_id(df), df['BindingDB Target Chain Sequence']):
            file.write(f">{id}|{name}\n{seq}\n")