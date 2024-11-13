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


def get_patent_info(patent_number):
    """
    Get information of a patent using google patents

    :patent_number: str, patent number

    return: list containing patent title, abstract, url, status, citations
    """
    url = f"https://patents.google.com/patent/{patent_number}/en"
    try:
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        title = soup.find("title").get_text()  # Find patent title
        abstract = soup.find("meta", {"name": "description"})[
            "content"
        ]  # Find patent abstract
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
