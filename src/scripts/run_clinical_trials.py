import os

import pandas as pd

from src.scripts import clinical_trials as cl

zinc_ids = pd.read_csv(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "data",
        "BindingDB_All.tsv",
    ),
    sep="\t",
    usecols=["ZINC ID of Ligand"],
).dropna()

zinc_ids_dict = {
    name: zinc_ids["ZINC ID of Ligand"].unique()[i::5]
    for i, name in zip(range(5), ["Wes", "Amelie", "Gregor", "Daphne", "Guillaume"])
}

cl.get_zinc_clinical_trial_data_for_all_ids(
    zinc_ids_dict["###### ECRIRE TON NOM A LA PLACE DE CA ######"]
)  # <-------- REMPLIR TON NOM ICI
