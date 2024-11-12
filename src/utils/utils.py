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
