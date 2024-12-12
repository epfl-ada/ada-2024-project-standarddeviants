import json
import os
from math import isnan
from typing import List, Union

import numpy as np
import pandas as pd
import requests
import seaborn as sns
from tqdm import tqdm


def get_zinc_clinical_trial_data(id: str) -> Union[List[Union[str, None]], str]:
    """Get ZINC clinical trial data from a ZINC ID

    Args:
        id (str): ZINC ID

    Returns:
        Union[List[Union[str, None]], str]: Clinical trials data
    """
    url = f"https://zinc.docking.org/substances/{id}/trials.json?count=all"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        return data
    except requests.RequestException as e:
        return np.nan


def get_zinc_clinical_trial_data_for_all_ids(
    ids: list[str], data_path: str = None
) -> dict[str, Union[List[Union[str, None]], str]]:
    """Get ZINC clinical trial data for a given list of ZINC IDs

    Args:
        ids (list[str]): ZINC IDs
        data_path (str, optional): The path with existing data or data to be created. Defaults to None.

    Returns:
        dict[str, Union[List[Union[str, None]], str]]: The clinical trials data
    """
    if data_path is None:
        data_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "data",
            "ZINC_references_trials.json",
        )

    zinc_data = {}
    if os.path.exists(data_path):
        with open(data_path, "r") as f:
            zinc_data = json.load(f)

    for i, id in enumerate(tqdm(ids)):
        if id not in zinc_data.keys():
            zinc_data.update({id: {"references": get_zinc_clinical_trial_data(id)}})
        if i != 0 and i % 100 == 0:
            with open(data_path, "w") as f:
                json.dump(zinc_data, f)
    return zinc_data


def load_clinical_trials_data(path_to_data: str = None) -> pd.DataFrame:
    """Load the clinicals trials data from the output of `get_zinc_clinical_trial_data_for_all_ids`.

    Args:
        path_to_data (str, optional): The path with existing data. Defaults to None.

    Returns:
        pd.DataFrame: The clinical trial data after some pre-processing
    """
    if path_to_data is None:
        path_to_data = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "data",
            "ZINC_references_trials.json",
        )
    df = pd.read_json(path_to_data).T
    df = df.rename(columns={"references": "trials"})
    df.index.name = "ZINC ID of Ligand"
    df_all = pd.DataFrame()
    for zinc_id, trials in df.iterrows():
        if not isinstance(trials["trials"], list) and isnan(trials["trials"]):
            continue
        trials_df = pd.DataFrame(trials["trials"])
        trials_df["ZINC ID of Ligand"] = zinc_id
        df_all = pd.concat([df_all, trials_df], axis=0)
    return df_all


def plot_clinical_counts(
    df: pd.DataFrame,
    ax,
    col: str,
    ylabel: str = None,
    title: str = None,
    top_n: int | None = None,
    yscale: str | None = None,
):
    """Custom coutplot for clinical trials data.

    Args:
        df (pd.DataFrame): Clinical trail dataframe - ouput from `load_clinical_trials_data`
        ax (matplotlib axis): matplotlib axis
        col (str): Column name in the df.
        ylabel (str, optional): The matplotlib ylabel. Defaults to None.
        title (str, optional): The matplotlib title. Defaults to None.
        top_n (int | None, optional): Shows the top_n first values in the column. Defaults to None.
        yscale (str | None, optional): The matplotlib yscale. Defaults to None.

    Returns:
        matplotlib axis: _description_
    """
    count_df = df[col].value_counts().reset_index()
    if top_n is not None:
        count_df.loc[top_n:, col] = "Other"
    sns.barplot(
        count_df,
        x=col,
        y="count",
        order=count_df[col][: top_n + 1] if top_n is not None else None,
        ax=ax,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if yscale is not None:
        ax.set_yscale(yscale)
    ax.set_xlabel(col.replace("_", " ").capitalize())
    return ax


if __name__ == "__main__":
    df = pd.read_csv(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "..",
            "data",
            "BindingDB_All.tsv",
        ),
        sep="\t",
        usecols=["ZINC ID of Ligand"],
    )
    data_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "data",
        "ZINC_references_trials.json",
    )
    ids = df["ZINC ID of Ligand"].value_counts().index

    get_zinc_clinical_trial_data_for_all_ids(ids, data_path)
