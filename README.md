# A Multi-Scale Analysis of Drug Discovery
## Abstract
Drug discovery is a decade-long process which often fails before clinical application and approval. Understanding its mechanisms will help avoid this failure, and highlight the patterns that yield success. For this reason, we propose a **multi-scale analysis** of drug discovery.

This analysis opens with **global trends**, such as contributions to drug discovery by countries, institutions and people, taking into account patents, publications and time. This will uncover foundational factors and topics, which guide the behaviour of research bodies. 
This will then allow investigating a level lower, to understand the patterns within, with a particular focus on researched disease areas, their model **organisms** and outcomes.
Finally, we explore the **molecular trends** and chemical workings of drug discovery themselves. Here, we aim to identify overrepresented chemical structures and particular binding kinetics, and how these change when targeting different protein classes. 
Overall, by travelling from kilometers to nanometers _(Ant-man x Jules Vernes)_, we will tell the evolving story of drug discovery.

## Research Questions
To perform the multi-scale analysis, we focus on a few particular research questions per scale.
### Global Scale
- What level of contributions are shown by countries and institutions to the global BindingDB dataset, and does this change with time?
- How do authors work together? Can we see any "author networks"?
- Leading to a lower scale, which diseases are covered the most in the dataset?
### Organism Scale
- In relation to the global scale, which organisms are prefered for each disease area?
- Leading to a lower scale, which targets/ligand pairs are studied in which organisms? What are their kinetics?
### Molecular Scale
- In relation to both the global and organism scale, what is and ideal protein target (per disease, per organism, etc.)?
- What are the chemical properties and structures of potent ligands?
- How do ligands' structures change between families of target proteins? How does this influence binding kinetics?

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
- **Amélie**: Global Scale analysis (Countries + Diseases)
- **Daphné**: Organism Scale analysis (Binding Kinetics + clustering (?))
- **Gregor**: Molecular Scale analysis (Potent structures + Kinetics)
- **Guillaume**: Molecular Scale analysis (Structure clustering)
- **Wesley**: Global Scale analysis (Authors)

