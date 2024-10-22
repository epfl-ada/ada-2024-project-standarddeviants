import requests
import json
import pandas as pd
import os

def get_extra_pub_data(doi:str)->list|str:
    """
    Gets metadata (journal, publication date and publisher) about a DOI.

    Args:
        doi (str): valid DOI link

    Returns:
        list: [journal, date, publisher]
    """
    url = f"https://api.crossref.org/works/{doi}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        #print(json.dumps(data, indent=4, sort_keys=True))
        #collected metadata
        journal = data.get('message', {}).get('container-title', [[None]])[0]
        publisher = data.get('message', {}).get('publisher', [[None]])
        publication_date = data.get('message', {}).get('issued', {}).get('date-parts', [[None]])[0]
        
        if publication_date or journal:
            return list((journal, publication_date[0], publisher))
        else:
            return "Neither publication date nor journal found"
    else:
        return f"Error: {response.status_code}"
    

def add_extra_pub_data(df:pd.DataFrame, filename:str="bindingDB_with_metadata.csv")->None:
    """
    Saves metadata of each unique DOI in a dataframe as a CSV file.
    I.e. the saved CSV file has n rows (excluding header), where n is 
    the number of unique DOIs of the dataframe. 
    
    Note: the CSV will generally be used to merge to the initial dataframe, on "Article DOI".

    Args:
        df (pd.DataFrame): dataframe containing at least a column specifically named "Article DOI", which contains DOI links to articles.
        filename (str): name of the file to save CSV to.
    """
    # group all identical DOIs together and get their pub year
    new_metadata = pd.DataFrame(df["Article DOI"].unique()).rename({0:"Article DOI"},axis=1)
    new_metadata[["journal", "year", "publisher"]] = new_metadata["Article DOI"].apply(lambda x: pd.Series(get_extra_pub_data(x)))

    # verify filename has .csv extension
    root, ext = os.path.splitext(filename)
    if not ext==".csv":
        ext='.csv'
        filename = root + ext

    # paths for save    
    base_dir = os.path.join(os.path.dirname(os.path.dirname( __file__ )))
    deriv_path = os.path.join(base_dir,"Derivatives")
    save_path = os.path.join(deriv_path,filename)

    # save in Derivatives folder
    if not os.path.exists(deriv_path):
        os.mkdir(deriv_path)
    new_metadata.to_csv(save_path)