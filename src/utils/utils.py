from thefuzz import fuzz, process
from tqdm import tqdm
from typing import Optional, Dict, List 
import pandas as pd


def lower(s):
    return str(s).lower()


# Helper function to find the best match
def fuzzy_merge(df1, df2, left_on, right_on, how, threshold=80):
    matches = []
    """
    Fuzzy merges two DataFrames by matching similar strings in specified columns.

    Parameters
    ----------
    df1 : pd.DataFrame
        The left DataFrame to merge.
    df2 : pd.DataFrame
        The right DataFrame to merge.
    left_on : str
        Column name in `df1` to match on.
    right_on : str
        Column name in `df2` to match on.
    how : str
        Type of merge ('left', 'right', 'inner', or 'outer').
    threshold : int, optional
        Minimum similarity score (0-100) for a match; default is 80.

    Returns
    -------
    pd.DataFrame
        A merged DataFrame with an added `fuzzy_match_<left_on>` column in `df1`, containing the 
        best match in `df2` if it meets the threshold.

    """
    
    for name in tqdm(df1[left_on]):
        # Get the best match with a score over the threshold
        match, score, _ = process.extractOne(
            name, df2[right_on], scorer=fuzz.token_sort_ratio
        )
        if score >= threshold:
            matches.append(match)
        else:
            matches.append(None)
    df1[f"fuzzy_match_{left_on}"] = matches
    # Merge on the original and the matched column
    return df1.merge(
        df2, left_on=f"fuzzy_match_{left_on}", right_on=right_on, how="left"
    )


def group_categories(
     name: str,
    in_mapping: Optional[dict[str, list[str]]],
    check_key_for_in_mapping: bool = True,
    endswith_mapping: Optional[dict[str, list[str]]] = None,
    priority_list: Optional[List[str]] = None
) -> str:
    """_summary_

    Args:
        nname (str): The name to be categorized.
        in_mapping (dict[str, list[str]] | None): Exact matches for categories. If None, no exact match is performed.
        check_key_for_in_mapping (bool, optional): Whether to check for exact matches in `in_mapping`. Default is True.
        endswith_mapping (dict[str, list[str]] | None, optional): Matches categories based on suffixes. Defaults to None.
        priority_list (List[str] | None, optional): Priority order for matching categories. Defaults to mapping order.

    Returns:
        str: _description_
    """
    name_lower = name.lower()

    # If a priority list is provided, use it to order the items in `in_mapping`
    if priority_list:
        ordered_mapping = [(k, in_mapping[k]) for k in priority_list if k in in_mapping]
    else:
        ordered_mapping = in_mapping.items()

    # Check in_mapping according to priority order
    for k, v in ordered_mapping:
        if check_key_for_in_mapping and (k.lower() in name_lower):
            return k
        for v_i in v:
            if v_i in name_lower:
                return k

    # Check endswith_mapping if no match was found in in_mapping
    if endswith_mapping is not None:
        for k, v in endswith_mapping.items():
            for v_i in v:
                if name_lower.endswith(v_i):
                    return k

    # Return original name if no category match is found
    return name


if __name__ == "__main__":
    in_mapping = {
        "Neurodegeneration": [
            "parkinson",
            "alzeihmer",
        ],
        "Cancer": [
            "leukemia",
            "oncogene",
        ],
        "Obesity": [],
        "Neoplasia": [],
        "Hirschsprung Disease": [],
        "Immunodeficiency": [],
        "Agammaglobulinemia": [],
        "Thrombocythemia": [],
        "Diabetes": [],
        "Cornelia de Lange syndrome": [],
        "QT syndrome": [],
    }

    endswith_mapping = {"Cancer": ["oma"]}

    f = lambda x: group_categories(
        x,
        in_mapping=in_mapping,
        check_key_for_in_mapping=True,
        endswith_mapping=endswith_mapping,
    )

    # df[... Class] = df[...].apply(f)


def count_classified_rows(df: pd.DataFrame) -> int:
    """
    Count how many rows in the DataFrame have been classified based on the `Target Name` 
    and 'Target Class' columns. A row is considered classified if the `Target Class` 
    is different from 'Target Name'.

    Args:
        df (pd.DataFrame): The dataframe containing the `Target Name` and `Target Class` columns.

    Returns:
        int: The count of classified rows.
    """
    if 'Target Name' not in df or 'Target Class' not in df:
        raise ValueError("DataFrame must contain 'Target Name' and 'Target Class' columns")
    
    # Count rows where 'Target Class' is not empty or NaN
    classified_count = ((df['Target Name'] != df['Target Class']) & (df['Target Class'].notna()))
    
    # Return the number of classified rows
    return classified_count.sum()


# print(count_classified_rows(...)) 