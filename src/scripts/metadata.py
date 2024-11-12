import json
import os
import pandas as pd
import requests
from typing import Union, List, Tuple 
from tqdm import tqdm 

def get_extra_pub_data(doi: str) -> Union[List[Union[str, None]], str]:
    """
    Gets metadata (journal, publication date, and publisher) about a DOI.

    Args:
        doi (str): Valid DOI link.

    Returns:
        Union[List[Union[str, None]], str]: A list containing [journal, date, publisher],
                                             or a string describing an error if the data
                                             cannot be retrieved.
    """
    url = f"https://api.crossref.org/works/{doi}"
    try:
        response = requests.get(url)
        response.raise_for_status()  
        data = response.json()
        
        # Safely extract the fields with default values if they are missing 
        journal = data.get("message", {}).get("container-title", [None])[0]
        publisher = data.get("message", {}).get("publisher", None)
        publication_date = data.get("message", {}).get("issued", {}).get("date-parts", [[None]])[0][0]

        # Return the collected data as a list
        if journal or publication_date or publisher:
            return [journal, publication_date, publisher]
        else:
            return "Neither publication date nor journal found"
    
    except requests.RequestException as e:
        # Handle network-related issues or invalid DOI requests
        return f"Error: {e}"
    except (IndexError, KeyError) as e:
        # Handle cases where expected fields are missing in the response structure
        return "Error in parsing response: missing fields in metadata"



def add_extra_pub_data(
    df: pd.DataFrame, filename: str = "bindingDB_with_metadata.csv", chunk_size: int = 1000
) -> None:
    """
    Saves metadata of each unique DOI in a dataframe as a CSV file in chunks.
    Each chunk will contain up to 5000 rows of DOI metadata.

    Args:
        df (pd.DataFrame): dataframe containing at least a column specifically named "Article DOI", which contains DOI links to articles.
        filename (str): name of the file to save CSV to.
        chunk_size (int): number of rows per saved chunk.
    """
    # Group all identical DOIs together and get their metadata
    unique_dois = pd.DataFrame(df["Article DOI"].unique(), columns=["Article DOI"])

    # Enable tqdm for pandas apply method
    tqdm.pandas(desc="Fetching metadata")

    # Verify filename has .csv extension
    root, ext = os.path.splitext(filename)
    if ext != ".csv":
        ext = ".csv"
        filename = root + ext

    # Define save path
    base_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)))
    deriv_path = os.path.join(base_dir, "data")
    save_path = os.path.join(deriv_path, filename)

    # Ensure "data" directory exists
    if not os.path.exists(deriv_path):
        os.mkdir(deriv_path)

    # Process data in chunks and save incrementally
    for start in range(0, unique_dois.shape[0], chunk_size):
        end = start + chunk_size
        chunk = unique_dois.iloc[start:end]

        # Fetch metadata for each DOI in the chunk
        chunk[["journal", "year", "publisher"]] = chunk["Article DOI"].progress_apply(
            lambda x: pd.Series(get_extra_pub_data(x))
        )

        # Append chunk to CSV file
        mode = 'w' if start == 0 else 'a' 
        header = (start == 0)
        chunk.to_csv(save_path, mode=mode, header=header, index=False)


if __name__ == "__main__":
    usecols = ["Article DOI"]
    df = pd.read_csv(
        r"/Users/.../code/ada-2024-project-standarddeviants/data/BindingDB_All.tsv",
        sep="\t",
        usecols=usecols,
    )
    add_extra_pub_data(df)
