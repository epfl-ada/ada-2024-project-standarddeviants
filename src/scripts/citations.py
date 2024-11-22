import re

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup


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


def calculate_H_index(disease_citations_df, citations_dict):
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
