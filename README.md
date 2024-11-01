# Breaking Good: How to make a successful drug
## Abstract
Drug discovery is a decade-long process which often fails before clinical application or approval. Understanding its mechanisms will help avoid this failure and lean towards success. For this reason, we propose an analysis on the biological, chemical and institutional indicators for success in drug discovery.

In particular, we will first show the different usage of model organisms to study ligands and targets. Thereafter, focusing on particular organisms, we are interested in discovering the underlying factors that may influence the outcome of drug discovery. In particular, we will investigate the biochemical and molecular features of ligands and targets, and how these change between diseases and target classes. Finally, we will link these findings to key indicators of success, such as publication citations, patents and drug status, to identify the best-in-class drugs and their traits that lead to success.


## Research Questions
#### Organisms
To guide this research, we will attempt to answer the following questions:
- Which target organisms are most frequently used in research
- How does distribution change depending on the the disease or type of protein targeted?

#### Protein Classes
- What are the most targeted protein classes in research and development? How does this relate to diseases?
- **What are the most "successful" target classes, in terms of publications, citations, patents and drug status?**

#### Molecular Features
- **What are the most "successful" ligands of different protein classes?**
- What are there molecular and biochemical features?
- Are these features shared accross multiple successful molecules?
  

## Datasets
In addition to the [BindingDB](https://www.bindingdb.org/rwd/bind/index.jsp) dataset, we will access the [DrugBank](https://en.wikipedia.org/wiki/DrugBank) dataset, to obtain information about disease areas covered by ligands. We will access information about a ligand in this dataset through its DrugBank ID, in particular its description and indications. From the retrieved text, we will use naural language processing tools to extract a list of diseases covered by each ligand.

To be continued...

## Methods
To be continued...

## Timeline and Team Organisation
### Milestones and Timeline
To be discussed...

### Team
Very loose proposal of roles (admin tasks still need assigning):
- TBD



## Conda Environment Setup

```shell
# Navigate to directory
cd .../ada-2024-project-standarddeviants

# Create a new conda environment with Python 3.11
conda create -n ada python=3.11.9

# Activate the environment
conda activate ada

# Install the required dependencies
pip install -r requirements.txt

# Set up pre-commit hooks
pre-commit install
```
