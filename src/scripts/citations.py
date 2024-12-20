import re

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
import json
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

def get_citations(doi):
    """
    Get number of citations based on a DOI

    :doi: str, DOI of article

    return: number of citations (int)
    """
    url = f"https://api.crossref.org/works/{doi}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        citations = data.get("message", {}).get("is-referenced-by-count", {})
        return citations
    else:
        print(f"Error: {response.status_code}")
        return None


def get_citations_per_disease(df, diseases_df, citations_dict):
    """Creates a dataframe containing articles and citations information for each disease class.

    Args:
        df (DataFrame): BindingDB dataframe
        diseases_df (DataFrame): df containing the Uniprot id under "UniProt (SwissProt) Primary ID of Target Chain" and the list of correponding disease classes under "Disease Classes".
        citations_dict (Dictionnary): for each article doi "doi" has the corresponding citation count "citation".

    Returns:
        DataFrame: contains the list of corresponding article DOIs, article count, total citations count, for each disease class.
    """
    disease_class_dict = {
        disease_class: set(
            diseases_df[
                diseases_df["Disease Classes"].apply(lambda x: disease_class in x)
            ]["UniProt (SwissProt) Primary ID of Target Chain"]
        )
        for disease_class in set(sum(diseases_df["Disease Classes"], []))
    }
    disease_class_df = pd.DataFrame(
        [
            (disease_class, uniprot_id)
            for disease_class, uniprot_ids in disease_class_dict.items()
            for uniprot_id in uniprot_ids
        ],
        columns=["Disease Classes", "UniProt (SwissProt) Primary ID of Target Chain"],
    )
    bindingDB = df.dropna(subset=["Article DOI"])
    merged_data = bindingDB.merge(
        disease_class_df,
        on="UniProt (SwissProt) Primary ID of Target Chain",
        how="inner",
    )
    disease_dois = (
        merged_data.groupby("Disease Classes")["Article DOI"].apply(set).reset_index()
    )
    disease_dois["Article Count"] = disease_dois["Article DOI"].apply(len)
    disease_dois["Total Citations"] = disease_dois["Article DOI"].apply(
        lambda dois: sum(citations_dict.get(doi, 0) or 0 for doi in dois)
    )
    disease_dois["Article DOI"] = disease_dois["Article DOI"].apply(list)
    disease_citations_df = disease_dois.sort_values(
        by="Article Count", ascending=False
    ).reset_index(drop=True)
    return disease_citations_df

def get_citations_per_target(names_df, citations_dict, mapped_names):
    """Creates a dataframe containing articles and citations information for each target class.

    Args:
        names_df (pd.DataFrame): BindingDB dataframe containing the "Article DOI",  "UniProt (SwissProt) Recommended Name of Target Chain", and "UniProt (TrEMBL) Submitted Name of Target Chain" columns.
        citations_dict (pd.DataFrame): for each article doi "doi" has the corresponding citation count "citation".
        mapped_names (pd.DataFrame): returned when using fonction targets.get_target_class(names_df=names_df).

    Returns:
        pd.DataFrame: contains the list of corresponding article DOIs and citations, for each target class.
    """
    merged = names_df.merge(mapped_names, left_index=True, right_index=True)
    merged = merged.dropna(subset="Article DOI")
    target_dois = (
        merged.groupby("UniProt (SwissProt) Recommended Name of Target Chain_y")["Article DOI"].apply(set).reset_index()
    )
    target_citations_exp = []
    for index, row in target_dois.iterrows():
        target_class = row['UniProt (SwissProt) Recommended Name of Target Chain_y']
        dois = row["Article DOI"]
        citations_list = []
        for doi in dois: 
            citation_count = citations_dict.get(doi, 0)  # 0 if DOI not found
            citations_list.append({"DOI": doi, "Citations": citation_count})
        target_citations_exp.append(
            {"Target Classes": target_class, "Citations": citations_list}
        )
    target_citations_df = pd.DataFrame(target_citations_exp)
    return target_citations_df

def calculate_H_index_diseases(disease_citations_df, citations_dict):
    """Calculates the H index for each disease class based on number of articles and citations.

    Args:
        disease_citations_df (DataFrame): should contain the list of corresponding article DOIs as "Article DOI" for each disease class ("Disease Classes").
        citations_dict (Dictionnary): for each article doi "doi" has the corresponding citation count "citation".

    Returns:
        DataFrame: contains calculated H-index for each disease class.
    """
    disease_citations_exp = []
    for index, row in disease_citations_df.iterrows():
        disease_class = row["Disease Classes"]
        dois = row["Article DOI"]
        citations_list = []
        for doi in dois:
            citation_count = citations_dict.get(doi, 0)  # 0 if DOI not found
            citations_list.append({"DOI": doi, "Citations": citation_count})
        disease_citations_exp.append(
            {"Disease Classes": disease_class, "Citations": citations_list}
        )
    disease_citations_exp_df = pd.DataFrame(disease_citations_exp)
    disease_h_index_exp = []
    for index, row in disease_citations_exp_df.iterrows():
        disease_class = row["Disease Classes"]
        citations_counts = [
            entry["Citations"]
            for entry in row["Citations"]
            if entry["Citations"] is not None
        ]
        if citations_counts:
            h_index = sum(
                x >= i + 1 for i, x in enumerate(sorted(citations_counts, reverse=True))
            )
        else:
            h_index = 0
        disease_h_index_exp.append(
            {"Disease Classes": disease_class, "H-Index": h_index}
        )
        disease_h_index_df = pd.DataFrame(disease_h_index_exp)
        disease_h_index_df = disease_h_index_df.sort_values(
            by="H-Index", ascending=False
        )
    return disease_h_index_df

def calculate_H_index_targets(target_citations_df):
    """Calculates the H index for each target class based on number of articles and citations.

    Args:
        target_citations_df (DataFrame): should contain the list of article DOIs and corresponding citations ("Citations") for each target class ("Target Classes").
    Returns:
        DataFrame: contains calculated H-index for each disease class.
    """
    target_h_index_exp = []
    for index, row in target_citations_df.iterrows():
        target_class = row["Target Classes"]
        citations_counts = [
            entry["Citations"]
            for entry in row["Citations"]
            if entry["Citations"] is not None
        ]
        if citations_counts:
            h_index = sum(
                x >= i + 1 for i, x in enumerate(sorted(citations_counts, reverse=True))
            )
        else:
            h_index = 0
        target_h_index_exp.append(
            {"Target Classes": target_class, "H-Index": h_index}
        )
        target_h_index_df = pd.DataFrame(target_h_index_exp)
        target_h_index_df = target_h_index_df.sort_values(
                by="H-Index", ascending=False
            )
    target_h_index_df
    return target_h_index_df

def calculate_patent_H_index_diseases(disease_citations_df):
    """Calculates the H index for each disease class based on number of patents citations.

    Args:
        disease_citations_df (DataFrame): should contain the list of patent number and corresponding citations ("Patent Citations") for each disease class ("Disease Classes").
    Returns:
        DataFrame: contains calculated H-index for each disease class.
    """
    disease_h_index_exp = [] 
    for index, row in disease_citations_df.iterrows():
        disease_class = row["Disease Classes"]
        citations_counts = [
            entry["Citations"]
            for entry in row["Patent Citations"]
            if entry["Citations"] is not None
        ]
        if citations_counts:
            h_index = sum(
                x >= i + 1 for i, x in enumerate(sorted(citations_counts, reverse=True))
            )
        else:
            h_index = 0
        disease_h_index_exp.append(
            {"Disease Classes": disease_class, "H-Index": h_index}
        )
        disease_h_index_df = pd.DataFrame(disease_h_index_exp)
        disease_h_index_df = disease_h_index_df.sort_values(
                by="H-Index", ascending=False
            )
    return disease_h_index_df

def calculate_patent_H_index_targets(targets_citations_df):
    """Calculates the H index for each target class based on number of patents citations.

    Args:
        targets_citations_df (DataFrame): should contain the list of patent number and corresponding citations ("Patent Citations") for each target class ("Target Classes").
    Returns:
        DataFrame: contains calculated H-index for each target class.
    """
    target_h_index_exp = []
    for index, row in targets_citations_df.iterrows():
        target_class = row["Target Classes"]
        citations_counts = [
            entry["Citations"]
            for entry in row["Citations"]
            if entry["Citations"] is not None
        ]
        if citations_counts:
            h_index = sum(
                x >= i + 1 for i, x in enumerate(sorted(citations_counts, reverse=True))
            )
        else:
            h_index = 0
        target_h_index_exp.append(
            {"Target Classes": target_class, "H-Index": h_index}
        )
        target_h_index_df = pd.DataFrame(target_h_index_exp)
        target_h_index_df = target_h_index_df.sort_values(
                by="H-Index", ascending=False
            )
    return target_h_index_df

def get_citations_per_year_diseases(disease_citations_df, citations_dict):
    """Get the total number of citations per year for each disease class.

    Args:
        disease_citations_df (DataFrame): should contain the list of corresponding article DOIs as "Article DOI" for each disease class ("Disease Classes").
        citations_dict (Dictionnary): for each article doi "doi" has the corresponding citation count "citation".

    Returns:
        DataFrame: contains the total number of citations per year for each disease class.
    """
    doi_metadata = pd.read_csv("src/data/metadata.csv").dropna()
    doi_year_dict = dict(zip(doi_metadata["Article DOI"], doi_metadata["year"]))
    disease_years = []
    for index, row in disease_citations_df.iterrows():
        disease_class = row["Disease Classes"]
        dois = row["Article DOI"]
        years = []
        citations_count = []
        for doi in dois:
            year = doi_year_dict.get(doi)
            citation_count = citations_dict.get(doi, 0)
            if year is not None:
                years.append(year)
                citations_count.append(citation_count)
            disease_years.append(
                {
                    "Disease Classes": disease_class,
                    "Publication Years": years,
                    "Citations": citations_count,
                }
            )
    disease_years_df = pd.DataFrame(disease_years)
    disease_years_df = disease_years_df.drop_duplicates(subset="Disease Classes")
    diseases_top_10 = disease_years_df.head(10)
    disease_years_expanded = []
    for _, row in diseases_top_10.iterrows():
        disease_class = row["Disease Classes"]
        years = row["Publication Years"]
        citations = row["Citations"]
        for year, citation in zip(years, citations):
            disease_years_expanded.append(
                {"Year": year, "Disease Classes": disease_class, "Citations": citation}
            )
    expanded_df = pd.DataFrame(disease_years_expanded)
    aggregated_df = expanded_df.groupby(["Year", "Disease Classes"], as_index=False)[
        "Citations"
    ].sum()
    return aggregated_df

def get_citations_per_year_targets(citations_list, doi_year_dict):
    """Get the total number of citations per year for each target class.

    Args:
        citations_list (list): for each article doi "doi" has the corresponding citation count "citation".
        doi_year_dict (dict): for each article doi has the corresponding year of publication.

    Returns:
        list: contains all the publication years corresponding to the input DOIs.
    """
    years = []
    for citation in citations_list:
        doi = citation.get('DOI')
        if doi in doi_year_dict:
            years.append(doi_year_dict[doi])
        else:
            years.append(None)
    return years

def timeseries_citations_diseases(citations_dict:dict, disease_dois:pd.DataFrame):
    """Creates a dataframe containing all the citations and year of publication for each article of a disease class.

    Args:
        citations_dict (dict): contains the article DOI and corresponding number of citations.
        disease_dois (pd.DataFrame): contains the disease class and the list of corresponding article DOIs.

    Returns:
        pd.DataFrame: contains all the citations and year of publication for each article of a disease class.
    """
    doi_metadata = pd.read_csv("../src/data/metadata.csv").dropna()
    doi_year_dict = dict(zip(doi_metadata["Article DOI"], doi_metadata["year"]))
    disease_years = []
    for index, row in disease_dois.iterrows():
        disease_class = row["Disease Classes"]
        dois = row["Article DOI"]
        years = []
        citations_count = []
        for doi in dois:
            year = doi_year_dict.get(doi)
            citation_count = citations_dict.get(doi, 0)
            if year is not None:
                years.append(year)
                citations_count.append(citation_count)
            disease_years.append(
                {
                    "Disease Classes": disease_class,
                    "Publication Years": years,
                    "Citations": citations_count,
                }
            )
    disease_years_df = pd.DataFrame(disease_years)
    disease_years_df = disease_years_df.drop_duplicates(subset="Disease Classes")
    return disease_years_df

def get_patent_info(patent_number):
    """
    Get information of a patent using google patents

    Args:
    patent_number: str, patent number

    Returns: list containing patent title, abstract, url, status, citations
    """
    url = f"https://patents.google.com/patent/{patent_number}/en"
    try:
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        title_ = soup.find("title")
        title = title_.get_text() if title_ else "Unknown Title"
        abstract_ = soup.find("meta", {"name": "description"})
        abstract = (
            abstract_["content"]
            if abstract_ and "content" in abstract_.attrs
            else "No abstract available"
        )
        status = soup.find("span", {"itemprop": "ifiStatus"})  # Find patent status
        status = status.get_text() if status else "N/A"
        families_citing_header = soup.find(
            "h2", string=re.compile("Families Citing this family")
        )  # Find number of times it was cited by "family to family citations"
        if families_citing_header:
            families_citing_text = families_citing_header.get_text()
            families_citing_count = re.search(r"\((\d+)\)", families_citing_text)
            families_citing_count = (
                families_citing_count.group(1) if families_citing_count else "0"
            )
        else:
            families_citing_count = "0"
        cited_by_header = soup.find(
            "h2", string=re.compile("Cited By")
        )  # Find number of times it was cited by others
        if cited_by_header:
            cited_by_text = cited_by_header.get_text()
            cited_by_count = re.search(r"\((\d+)\)", cited_by_text)
            cited_by_count = cited_by_count.group(1) if cited_by_count else "0"
        else:
            cited_by_count = "0"
        return {
            "title": title,
            "abstract": abstract,
            "url": url,
            "status": status,
            "families citing": families_citing_count,
            "cited by": cited_by_count,
        }
    except requests.RequestException as e:
        print(f"Error fetching patent details: {e}")
        return None

def load_patents():
    """Loads the patents.json file and creates a pd.DataFrame with the relevant columns easy to use.

    Returns:
        pd.DataFrame: contains all the useful patent information.
    """
    with open("src/data/patents.json", "r") as f:
        patents = json.load(f)
    patents_df = pd.DataFrame(
        [
            {
                "patent": patent["patent"],
                "title": patent["info"].get("title", np.nan),
                "abstract": (patent["info"].get("abstract", "") or "").strip(),
                "url": patent["info"].get("url", np.nan),
                "status": patent["info"].get("status", np.nan),
                "families citing": int(patent["info"].get("families citing", 0) or 0),
                "cited by": int(patent["info"].get("cited by", 0) or 0),
            }
            for patent in patents
            if isinstance(patent, dict) and isinstance(patent.get("info"), dict)
        ]
    )
    return patents_df

def get_patent_citations_per_disease(df, diseases_df, patents_df):
    """Calculate total number of patent citations for each corresponding disease class.

    Args:
        df (pd.DataFrame): BindingDB dataframe containing the columns "Patent Number" and "UniProt (SwissProt) Primary ID of Target Chain".
        diseases_df (pd.DataFrame): should contain disease classes and associated Uniprot IDs.
        patents_df (pd.DataFrame): should contain patents and their citation counts.

    Returns:
        pd.DataFrame: with disease classes and their corresponding patents, total patent citations, and patent information, sorted by total citation count.
    """
    disease_class_dict = {
        disease_class: set(
            diseases_df[
                diseases_df["Disease Classes"].apply(lambda x: disease_class in x)
            ]["UniProt (SwissProt) Primary ID of Target Chain"]
        )
        for disease_class in set(sum(diseases_df["Disease Classes"], []))
    }
    disease_class_df = pd.DataFrame(
        [
            (disease_class, uniprot_id)
            for disease_class, uniprot_ids in disease_class_dict.items()
            for uniprot_id in uniprot_ids
        ],
        columns=["Disease Classes", "UniProt (SwissProt) Primary ID of Target Chain"],
    )
    bindingDB = df.dropna(subset=["Patent Number"])
    merged_data = bindingDB.merge(
        disease_class_df,
        on="UniProt (SwissProt) Primary ID of Target Chain",
        how="inner",
    )
    disease_patents = (
        merged_data.groupby("Disease Classes")["Patent Number"].apply(set).reset_index()
    )
    disease_patents["Total Citations"] = disease_patents["Patent Number"].apply(
        lambda patents: sum(
            patents_df.loc[patents_df["patent"].isin(patents), "total citations"]
        )
    )
    def get_patent_citations(patents): 
        citations = []
        for patent in patents:
            citation_data = patents_df.loc[patents_df["patent"] == patent, "total citations"]
            citation_count = citation_data.iloc[0] if not citation_data.empty else 0
            citations.append({"Patent Number": patent, "Citations": citation_count})
        return citations
    disease_patents["Patent Citations"] = disease_patents["Patent Number"].apply(get_patent_citations)
    disease_citations_df = disease_patents.sort_values(
        by="Total Citations", ascending=False
    ).reset_index(drop=True)
    return disease_citations_df


def get_patent_citations_per_targets(df, mapped_names, patents_df):
    """Gets the patent citations for each target class.

    Args:
        df (pd.DataFrame): BindingDB dataframe.
        mapped_names (pd.DataFrame): returned when using fonction targets.get_target_class(names_df=names_df).
        patents_df (pd.DataFrame): for eahc "Patent Number" should contain the corresponding "total citations".

    Returns:
        pd.DataFrame: contains the columns "Target Classes" and "Citations" with a list of patents and corresponding citations for each target class.
    """
    merged = df.merge(mapped_names, left_index=True, right_index=True)
    merged = merged.dropna(subset="Patent Number")
    target_dois = (
        merged.groupby("UniProt (SwissProt) Recommended Name of Target Chain_y")["Patent Number"].apply(set).reset_index()
    )
    patents_dict = patents_df.set_index('patent')['total citations'].to_dict()
    target_citations_exp = []
    for index, row in target_dois.iterrows():
        target_class = row['UniProt (SwissProt) Recommended Name of Target Chain_y']
        patents = row["Patent Number"]
        citations_list = []
        for patent in patents: 
            citation_count = patents_dict.get(patent, 0)  # 0 if DOI not found
            citations_list.append({"Patent": patent, "Citations": citation_count})
        target_citations_exp.append(
            {"Target Classes": target_class, "Citations": citations_list}
        )
    target_patent_citations = pd.DataFrame(target_citations_exp)
    return target_patent_citations

def plot_patents_disease(disease_patents_df, patents_h_index):
    """Creates three subplots, one for patent count per disease class, another for patent citations per disease class, and a third one for h-index metric per disease class.

    Args:
        disease_patents_df (pd.DataFrame): should contain for each disease class, the total patent citations and patents count.
        patents_h_index (pd.DataFrame): should contain for each disease class, the patent h-index value previously calculated.
    """
    palette_left = sns.color_palette("dark", len(disease_patents_df[0:10]))
    disease_classes = disease_patents_df['Disease Classes'][0:10].values
    color_dict = {disease_classes[i]: palette_left[i] for i in range(len(disease_classes))}
    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(2, 3, width_ratios=[2, 2, 4], height_ratios=[1, 1], figure=fig)

    ax_top = fig.add_subplot(gs[0, 0:2])
    ax_bottom = fig.add_subplot(gs[1, 0:2], sharex=ax_top)
    ax_hindex = fig.add_subplot(gs[:, 2])

    # 1st plot: 
    sns.barplot(
        data=disease_patents_df[0:10],
        x="Disease Classes",
        y="Patents Count",
        palette="dark",
        ax=ax_top,
    )
    ax_top.set_title("Distribution of Patents Count per Disease Class")
    ax_top.set_ylabel("Patents Count")
    ax_top.set_xlabel('')
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    # 2nd plot: 
    sns.barplot(
        data=disease_patents_df[0:10],
        x="Disease Classes",
        y="Total Citations",
        palette="dark",
        ax=ax_bottom,
    )
    ax_bottom.set_title("Distribution of Total Citations per Disease Class")
    ax_bottom.set_ylabel("Total Citations")
    ax_bottom.set_xlabel("Disease Classes (Top 10)")
    ax_bottom.set_xticklabels(ax_bottom.get_xticklabels(), rotation=90)

    # 3rd plot:
    # Adjust the color palette for the H-Index plot to match the disease classes
    new_palette = []
    new_disease_classes = patents_h_index['Disease Classes'].head(10).values
    used_colors = set(color_dict.values())
    for disease_class in new_disease_classes:
        if disease_class in color_dict:
            new_palette.append(color_dict[disease_class])
        else:
            new_color = sns.color_palette("Set2", len(new_disease_classes))[list(new_disease_classes).index(disease_class)]
            new_palette.append(new_color)

    sns.barplot(
        x="Disease Classes",
        y="H-Index",
        data=patents_h_index.head(10),
        palette=new_palette,
        ax=ax_hindex,
    )
    ax_hindex.set_title("Patent H-Index per Disease Class (Top 10)")
    ax_hindex.set_xlabel("Disease Classes")
    ax_hindex.set_ylabel("H-Index")
    ax_hindex.set_xticklabels(ax_hindex.get_xticklabels(), rotation=90)
    plt.tight_layout()
    plt.show()


def aggregate_citations_by_year_target(target_citations_df):
    """Aggregates citations by year and target class.

    Args:
        target_citations_df (pd.DataFrame): contains the column "Citations" with a list of citation counts (for each DOI),
                                                    the column "Publication Years" with a list of publication years (corresponding to each citation in Citations),
                                                    and the column "Target Classes" with the corresponding target class.

    Returns:
        pd.DataFrame: dataframe with aggregated citation counts, containing the column "Year" with the publication year, "Target Classes" with the target class, and "Citations" with the total citations count for the corresponding year and target class. 
    """
    expanded_rows = []
    for index, row in target_citations_df.iterrows():
        citations = row["Citations"]
        years = row["Publication Years"]
        target_class = row["Target Classes"]
        for year, citation in zip(years, citations):
            expanded_rows.append({"Year": year, "Target Classes": target_class, "Citations": citation})
    expanded_df = pd.DataFrame(expanded_rows)
    def extract_citations(citation_dict):
        if isinstance(citation_dict, dict) and 'Citations' in citation_dict:
            return citation_dict['Citations']
        else:
            return 0
    expanded_df["Citations"] = expanded_df["Citations"].apply(extract_citations)
    expanded_df["Citations"] = pd.to_numeric(expanded_df["Citations"], errors='coerce')
    aggregated_df = expanded_df.groupby(["Year", "Target Classes"], as_index=False)["Citations"].sum()
    return aggregated_df