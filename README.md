# Mapping the Molecular Path to Health: *A causal analysis of success in drug discovery*

## Abstract
Drug discovery is a decade-long process which often fails before clinical application or approval. Understanding its mechanisms can help avoid this failure and lean towards success. For this reason, we propose an analysis on the molecular, chemical and disease-related indicators for success in drug discovery.

More particularly, we first investigate the biochemical and molecular features of ligands, to allow their classification into different subtypes of drugs. These are subsequently compared with respect to the corresponding targeted protein classes and disease types. Concurrently, we develop success metrics, derived from publication citations, patents, clinical trials and drug status. We further analyze these metrics over time, offering a more dynamic depiction of success. Ultimately, we aim to causally link the uncovered biological and chemical features to the success metrics, in order to contribute to a more comprehensive landscape of drug discovery.

## Research Questions
To guide the research, we propose answering the following questions, focused on developing success metrics derived from citations, patents, prescriptions, and clinical trials to complement BindingDB.

#### **Molecular Features and Drug Success**:
- Do specific biochemical features of ligands (binding kinetics, structures …) predict their success metrics?
- Do certain properties like molecular weight, hydrophobicity, or binding affinity have a strong predictive relationship with success metrics?
- Is there a difference in success between competitive and non-competitive ligands?
- Do biochemical rules of thumb (eg. Lipinski’s rule) correlate with the derived success metrics?

#### **Hypothesis on Protein Class Success**:
- Do certain protein classes have a causal link to better success rates?
- In relation to disease specific outcomes, how and why are certain protein classes (e.g. GPCRs or kinases) preferentially targeted ? What class-specific molecular features might be favoured?

#### **Disease-Specific Outcomes**:
- Is there a causal relationship between the targeted diseases and success indicators?

#### **Temporal Trends of Success Metrics**:
- Do the success metrics change over time? What external factors may influence these?
- Do periods of intensified research activity, in certain protein classes or diseases, lead to higher success rates?
####  **Institutional Influence on Drug Discovery Outcomes**:
- Does the geographical location and institutional affiliation of research groups have a causal impact on the success rate of drug candidates?

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
A key pillar of this project is defining success metrics and understanding their limitations. Here, we outline current and planned metrics.
#### **Current metrics**:
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

#### **Future metrics to be explored**:
- Patents
    - Indicates the importance of drugs developed outside of academia
    - Limited by external confounding factors, such as institution’s influence and visibility
- Binding affinities and kinetics (eg. IC50 and Ki)
    - Indicates biochemical success, based on previous biological and chemical knowledge of optimal binding kinetics for particular classes of drugs
    - Limited to certain classes and by current state of biochemical knowledge

### Analysis
- Step 1: Data exploration
- Step 2: Data enrichment
	- 2.1: Uniprot data retrieval on diseases.
	- 2.2: Crossref API data retrieval on years and citations.
	- 2.3: DrugBank data retrieval on ligands.
    - 2.4: Medicaid and NDC Directory data retrieval on prescriptions.
	- 2.5: Zinc data retrieval on clinical trials.
	- 2.6: Google patents data retrieval on patents.
    - 2.7: Regression models to impute missing data, especially for binding kinetics.
- Step 3: Analysis and addressing research questions
    - 3.1: Develop success metrics as described above and review literature to determine optimal binding kinetics per disease class of interest.
    - 3.2: Analyse biochemical features, by clustering high-dimensional molecular data through dimensionality reduction techniques.
    - 3.3: Investigate success metrics’ relationships with biological and chemical features. We propose regression analysis with propensity score matching to control confounders. Success metrics, influenced by external factors, will be modelled as dependent variables, with protein class, ligand features, disease class, and binding kinetics as independent variables to establish causality.
    - 3.4: Access to publication dates enables time series analysis, revealing the dynamics of ligand research.

- Step 4: Create the data story webpage.

## Timeline and Team Organisation

![Dessin sans titre (1)](https://github.com/user-attachments/assets/66bc1b11-d585-4c36-a902-7099d5d8df7f)


- White: All
- Blue: Amélie and Wesley
- Green: Guillaume and Gregor
- Purple: Daphné and Gregor
- Orange: Wesley

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
