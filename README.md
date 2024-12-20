# [Mapping the Molecular Path to Health: *An analysis of success in drug discovery*](https://epfl-ada.github.io/ada-2024-project-standarddeviants/)


## Abstract
Drug discovery is a decade-long process which often fails before clinical application or approval. Understanding its mechanisms can therefore help avoid this failure and lean towards success. For this reason, we propose an analysis on the molecular, chemical and disease-related indicators for success in drug discovery. In particular, we show how distinct structural and chemical features, such as higher molecular weights, higher maximum absolute partial charges, and lower number of aliphatic carbocycles increase the efficacy of drugs. In practice, this is achieved by branching out and exploring the structural space by large screenings of slight modifications to an initial academic publication. This then brings new patents, new features and potentially future clinical trials, although these are lengthy in time. At a higher level, although historically established targets, such as CDKs, RTKs and growth factor receptors are still highly influential, there exist emerging targets, like neurotransmitter and hormone receptors. Globally, this analysis maps the current state of drug discovery, focusing on health issues of the 21st century, such as cancer, immunodeficiency, inflammation and neurodegeneration.

## Full Analysis
The detailed analysis is available on [our website](https://epfl-ada.github.io/ada-2024-project-standarddeviants/), with further details, derivations and supplementary material [here](/results.ipynb).

## Research Questions
To guide the research above, we answer three distinct research questions, all focused on success metrics derived from citations, patents, prescriptions, and clinical trials. Sub-questions help in formulating answers to this research.

#### **What are the relationships between targeted diseases and success in drug discovery?**
- What diseases are most influential? Which have the largest success?
- How do these metrics change over time?
  
#### **What are the relationships between targeted protein classes and success in drug discovery?**
- Are specific target protein linked to better success metrics?
- How are certain protein classes linked to disease classes identified as highly influential?
  
#### **What molecular features determine success in drug discovery?**
- Do specific biochemical features of ligands (binding kinetics, structures …) predict their success metrics?
- Do certain properties like molecular weight, hydrophobicity, or binding affinity have a strong predictive relationship with success metrics?
- How does institutional influence affect the research for new drug structures?
- What is the research strategy linked to developing new succesful structures?


## Datasets
- [BindingDB](https://www.bindingdb.org/rwd/bind/index.jsp)
[DrugBank](https://en.wikipedia.org/wiki/DrugBank)
    - Accessed via a download.
    - Used to retrieve generic names of commercialised drugs. One of the steps to link prescription data to the dataset.
- [ZINC](https://zinc.docking.org/)
    - Accessed with the Zinc Docking API : `https://zinc.docking.org/substances/{id}/trials.json?count=all` where id is the ZINC ID of the ligand found within BindingDB
    - Retrieving all clinical trial data for all ZINC IDs would take ~120 hours and is planned for phase 3. In phase 2, we processed the first 3000 most common ZINC IDs to preview the clinical trial data.
- [UNIPROT](https://www.uniprot.org/uniprotkb)
    - Accessed with the uniprot api
    - Contains information about disease areas covered by ligands
- [Google Patents](https://patents.google.com/)
    - Accessed to gather information about patents
- [National Drug Code (NDC) Directory](https://www.fda.gov/drugs/drug-approvals-and-databases/national-drug-code-directory)
    - Accessed via a download.
    - Used to retrieve brand names of commercialised drugs by merging them to the generic names in DrugBan, linking prescription data to the dataset.
- [Medicaid State Drug Utilization Data](https://www.medicaid.gov/medicaid/prescription-drugs/state-drug-utilization-data/index.html)
    - Accessed via a download.
    - Used to retrieve prescription data to the dataset by merging it to the brand names in the NDC directory.
- Article DOI Metadata
    - Accessed with the [crossref API](https://api.crossref.org/swagger-ui/index.html) to access publication metadata, based on a DOI (year of publication, journal, publisher, number of citations, authors, …)
## Methods
### Metrics
Key results of this study are based on success metrics, each with their limitations. Here, we outline these metrics.
- Number of publications related to a ligand
- Number of citations of initial publication of a ligand
    - Indicates the relevance of research on different target-ligand interactions
    - Limited by external confounding factors, such as author renown and affiliation, number of targets studies per publication, disease area, etc.
- H-index when considering each target or disease class as an entity
    - Compares the overall productivity and impact across different classes
    - May be skewed by duration of research activity, favouring more established fields over emerging ones
- Number and phase of clinical trials associated with each ligand.
    - Indicates the advancement of a drug in clinical development
    - Does not account for the quality or success rate of the trial.
- Number of prescriptions of the drug.
    - We only use data provided by medicaid about drug utilisation in the US, creating bias.
- Number of patent publications
- Number of citations linked to a patent
    - Indicates the importance of drugs developed outside of academia
    - Limited by external confounding factors, such as institution’s influence and visibility
- Binding affinities and kinetics (notably IC50)
    - Indicates biochemical success, based on biological of optimal binding kinetics for particular classes of drugs
    - IC50 is the most widely used measure of a drug's efficacy, indicating the concentration of a drug needed to inhibit a process by half.
    - Limited to certain classes and by current state of biochemical knowledge

### Analysis
#### Data exploration
We first explore the data and determine which methods of analysis, data transformations, and cleaning is necessary. Thereafter, we prepare the main data for all subsequemt analyses.

#### Data enrichment
We enrich the main BindingDB dataset, by retrieving the data presented above.
- Uniprot data retrieval on diseases
- Crossref API data retrieval on years and citations
- DrugBank data retrieval on ligands
- Zinc data retrieval on clinical trials
- Google patents data retrieval on patents

#### Methods for addressing research questions
We first develop the success metrics presented above, and review biochemical litereature per disease class of interest. Thereafter, we apply the metrics and rank both disease classes and target protein. We thereafter investigate these rankings and interpret the findings, also influenced by change over time, obtained through the publication date of each publication. In particular, this allows to identify two target proteins of interest (RET and D2R), due to their large implications in important diseases. These targets are then used as case studies for molecular feature analyses.

To analyse molecular features, we embed the molecular structures given by the SMILES data in 2048 bit vectors. These fingerprints capture structural information about a molecule, such as the presence of particular atoms, bonds, and functional groups. Thereafter, we group similar structures by reducing the dimensionality of data by principal component analysis (PCA). The obtained data points are grouped by K-means clustering, to identify families of ligands per target. Thereafter, we overlay the previously developed metrics, such as publication and patent citations, clinical trials and ligand efficacy, and identify how these clusters change over time.

Finally, we extract molecular features and structural elements (number of bonds, rings, partial charges, etc.) using [RDKit](https://www.rdkit.org/docs/index.html). Using IC50 as a measure of efficacy, we implement two classification models (multinomial logistic regression and random forests) to identify the predictive power and importance of the different chemical features.


## Timeline and Team Organisation
<img width="1073" alt="timeline" src="https://github.com/user-attachments/assets/5e54a7d1-ade2-4202-ac67-712335d4c61b" />
- Blue: Amélie
- Purple: Guillaume
- Green: Guillaume and Gregor
- Orange: Wesley and Amélie
- White: Wesley and Daphné
- Grey: Gregor

## Conda Environment Setup

```shell
# Navigate to directory
cd .../ada-2024-project-standarddeviants

# Create a new conda environment with Python 3.11
conda create -n ada python=3.11.9

# Activate the environment
conda activate ada

# Install the required dependencies
pip install -r pip_requirements.txt

# Set up pre-commit hooks
pre-commit install

# Install local src package
pip install -e . --use-pep517
```

## Repository Structure

```shell
.
├── README.md
├── data
│   ├── README.md # Explains how to get BindingDB_All.tsv
│   └── BindingDB_All.tsv
├── pip_requirements.txt
├── .gitignore
├── .pre-commit-config.yaml # For formatting code with black on every push
├── results.ipynb
├── setup.py # Setup file in order to build the local src package with pip
├── src
│   ├── __init__.py
│   ├── data # Collected external data
│   │   ├── README.md # Explains how the data was generated
│   │   ├── UniprotID_disases.json
│   │   ├── ZINC_references_trials.json
│   │   ├── citations.json
│   │   ├── metadata.csv
│   │   └── prescription_per_drugbank_id.csv
│   ├── model # Empty
│   ├── scripts
│   │   ├── __init__.py
│   │   ├── citations.py
│   │   ├── clinical_trials.py
│   │   ├── data_description.py
│   │   ├── disease_plotting.py
│   │   ├── geography.py
│   │   ├── metadata.py
│   │   ├── smiles.py
│   │   ├── targets.py
│   │   └── uniprot.py
│   └── utils
│       ├── __init__.py
│       └── utils.py
└── tests # Various notebooks where preliminary EDA was done
    ├── ...
    └── ...
```
