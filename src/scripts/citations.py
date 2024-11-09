import numpy as np
import pandas as pd
import requests


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
