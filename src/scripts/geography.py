import re

import geopandas as gpd
import numpy as np
import pandas as pd
import pycountry
import requests


def get_location(doi):
    """
    Get affiliated institution location based on a DOI

    :doi: str, DOI of article

    return: location (str)
    """
    url = f"https://api.crossref.org/works/{doi}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        authors = data.get("message", {}).get("author", [])
        if authors:
            main_author = authors[-1]
            affiliation_info = main_author.get("affiliation", [])
            if affiliation_info:
                location = affiliation_info[0].get("name", "No location found")
                return location
            else:
                return "No affiliation found"
        else:
            return "No authors found"

    else:
        return f"Error: {response.status_code}"


def add_location(df):
    """
    Add affiliated institution location to dataframes

    :df: dataframe containing a column of DOIs TITLED "Article DOI"

    return: new Pandas.Dataframe with an additional "Affiliated Location" column
    """
    # group all identical DOIs together and get their affiliated institution location
    pub_counts = df["Article DOI"].value_counts().to_frame().reset_index()
    pub_counts["Affiliated Location"] = pub_counts["Article DOI"].apply(get_location)

    # final dataframe with date of publications
    return df.merge(pub_counts, on="Article DOI", how="left")


def extract_country(affiliation):
    """
    Extracts the country from an affiliation address

    :affiliation: str, containing the address of interest

    return: the country (str)
    """
    us_states = {
        "Alabama",
        "Alaska",
        "Arizona",
        "Arkansas",
        "California",
        "Colorado",
        "Connecticut",
        "Delaware",
        "Florida",
        "Georgia",
        "Hawaii",
        "Idaho",
        "Illinois",
        "Indiana",
        "Iowa",
        "Kansas",
        "Kentucky",
        "Louisiana",
        "Maine",
        "Maryland",
        "Massachusetts",
        "Michigan",
        "Minnesota",
        "Mississippi",
        "Missouri",
        "Montana",
        "Nebraska",
        "Nevada",
        "New Hampshire",
        "New Jersey",
        "New Mexico",
        "New York",
        "North Carolina",
        "North Dakota",
        "Ohio",
        "Oklahoma",
        "Oregon",
        "Pennsylvania",
        "Rhode Island",
        "South Carolina",
        "South Dakota",
        "Tennessee",
        "Texas",
        "Utah",
        "Vermont",
        "Virginia",
        "Washington",
        "West Virginia",
        "Wisconsin",
        "Wyoming",
    }
    if not isinstance(affiliation, str):
        return None
    cleaned_affiliation = re.sub(r"\b[A-Z]{1,2}\d{1,4}[A-Z]{0,2}\b", "", affiliation)
    for (
        state
    ) in (
        us_states
    ):  # checks if US state present and map it to USA (some addresses from USA only have the state and not the country)
        if state in cleaned_affiliation:
            return "United States of America"
    # Matches countries in pyountry database (standardized country names in english):
    countries = []
    for country in pycountry.countries:
        if country.name in cleaned_affiliation or (
            country.alpha_2 in cleaned_affiliation
            or country.alpha_3 in cleaned_affiliation
        ):
            countries.append(country.name)
    if countries:  # returns the first matching country or None if no match
        return countries[0]
    return None
