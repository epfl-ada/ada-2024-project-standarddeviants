from thefuzz import fuzz, process
from tqdm import tqdm


def lower(s):
    return str(s).lower()


# Helper function to find the best match
def fuzzy_merge(df1, df2, left_on, right_on, how, threshold=80):
    matches = []
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
    in_mapping: dict[str, list[str]] | None,
    check_key_for_in_mapping: bool = True,
    endswith_mapping: dict[str, list[str]] | None = None,
) -> str:
    """_summary_

    Args:
        name (str): _description_
        in_mapping (dict[str, list[str]] | None): _description_
        check_key_for_in_mapping (bool, optional): _description_. Defaults to TRue.
        endswith_mapping (dict[str, list[str]] | None, optional): _description_. Defaults to None.

    Returns:
        str: _description_
    """
    name_lower = name.lower()

    for k, v in in_mapping.items():
        if check_key_for_in_mapping and (k.lower() in name_lower):
            return k
        for v_i in v:
            if v_i in name_lower:
                return k

    if endswith_mapping is not None:
        for k, v in endswith_mapping.items():
            for v_i in v:
                if name_lower.endswith(v_i):
                    return k

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
