import re

import geopandas as gpd
import numpy as np
import pandas as pd
import pycountry
import requests


def get_location(doi):
    """
    Retrieves the location of the affiliated institution for the main author of a given article based on its DOI.

    Parameters:
        doi (str): The DOI of the article.

    Returns:
        str: The name of the main author's affiliated institution or an error message if not found.
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
    Adds a column with affiliated institution locations to a DataFrame containing article DOIs.

    Parameters:
        df (pd.DataFrame): DataFrame with a column "Article DOI" containing DOIs of articles.

    Returns:
        pd.DataFrame: DataFrame with an added "Affiliated Location" column.
    """
    # group all identical DOIs together and get their affiliated institution location
    pub_counts = df["Article DOI"].value_counts().to_frame().reset_index()
    pub_counts["Affiliated Location"] = pub_counts["Article DOI"].apply(get_location)

    # final dataframe with date of publications
    return df.merge(pub_counts, on="Article DOI", how="left")


def extract_country(affiliation):
    """
    Extracts the country from an institution's affiliation address, including mappings for US states.

    Parameters:
        affiliation (str): The affiliation address containing location information.

    Returns:
        str or None: The country name, "United States of America" if a US state is detected,
                     or None if no country could be determined.
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
