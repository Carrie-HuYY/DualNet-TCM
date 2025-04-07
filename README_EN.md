# DualNet-TCM(Traditional Chinese Medicine Integrative Pharmacology Database)

***DualNet-TCM***(TCM (Traditional Chinese Medicine) + Dual (Dual) + Net (Network Structure) Integrative Pharmacology):
From the perspective of network pharmacology and traditional Chinese medicine theory, construct a bidirectional network of "Syndrome Differentiation - Formula - TCM - Components - Targets."

Simultaneously, develop diverse technologies: 1. Systematically screen and classify candidate proteins or genes to precisely narrow down the target range requiring experimental validation. 2. Incorporate safety assessment screening algorithms to assist in constructing a TCM toxicology network. 3. Text-based syndrome differentiation prediction to aid in drug repurposing and target discovery.

- [1. Introduction](#Introduction)

![figure1](#README_pictures/figure1.tiff)

- [2. Installation](#Installation)

- [3. Usage](#Usage)
  - [3.1 from_SD：Starting from Syndrome Differentiation](#from_sd)
  - [3.2 from_TCM_or_Formula：Starting from TCM/Formula](#from_tcm_or_formula)
  - [3.3 from_proteins：Starting from Targets](#from_protein)

- [4. Results Display](#Results Display)
  - [4.1 Interactive Graphs Based on pyecharts](#out_graph)
  - [4.2 Files for Cytoscape Visualization](#out_for_cytoscape)
  - [4.3 Excel Sheets for Saving Results](#out_for_excel)
  - [4.4 Query Results for Target Research Progress](#research_status_test)
  - [4.5 Safety Assessment Results](#safety_research)

- [5. Notes](#Notes)
  - [5.1 Regarding Python Version](#1Regarding Python Version)
  - [5.2 Regarding Data Download](#2Regarding Data Download)

## Introduction



## Installation

You can install DualNet-TCM using pip:

```cmd
pip install DualNet-TCM
```

Requirements include (Python version needs to be ≥3.9, others are not mandatory):

```cmd
pandas~=1.3.5
rdkit~=2023.3.2
tqdm~=4.67.1
pyecharts~=2.0.7
numpy~=1.21.6
```

## Install

### `from_SD`
Start from TCM syndrome differentiation to construct TCM integrative pharmacology.

```python
from DualNet-TCM import analysis
analysis.from_SD(
  sd, 
  score=990, 
  out_graph=False, 
  out_for_cytoscape=False, 
  out_for_excel=True,
  research_status_test=False,
  safety_research=False,
  re=True, 
  path='results/'
)
```

`from_SD`requires a mandatory parameter sd：any combinable data type that can use in to determine if an element is in the syndrome differentiation dataset, storing the IDs of the syndrome differentiations to be analyzed，
storing the IDs of the syndrome differentiations to be analyzed, such as ['DNS001'], Specific ID numbers and SD information can be found in the[SD Syndrome Differentiation Dataset].(/DualNet-TCM/Data/SD.xlsx)

`from_SD`Optional parameters：
- `score`: int type，in the [Chemical_Protein_Links Dataset](/DualNet-TCM/Data/Chemical_Protein_Links.xlsx), 
only records with a combined_score greater than or equal to scorewill be filtered out, default is 990；
- `out_graph`: boolean type, whether to output an ECharts-based interactive network visualization graph in HTML format, default is `False`；
- `out_for_cytoscape`: boolean type, whether to output files for Cytoscape visualization, default is `False`；
- `out_for_excel`: boolean type, whether to output results to an Excel file, default is `True`；
- `research_status_test`:  boolean type, whether to evaluate the research value of the obtained targets, default is `False`；
- `safety_research`: boolean type, whether to conduct safety research on the obtained targets, default is `False`；
- `re`: boolean type, whether to return the original analysis results (syndrome differentiation, formula, TCM, compounds (TCM components), proteins (targets), and their connection information), default is True。If re is True，the function will return the results sd,sd_formula_links.xlsx, formula,formula_tcm_links,tcm,tcm_chem_links,chem、
chem_protein_links and proteins, all of which arepd.DataFrame types, storing syndrome differentiation information, syndrome differentiation-formula connection information, formula information, formula-TCM connection information, TCM information, TCM-compound (TCM component) connection information, compound (TCM component) information, compound (TCM component)-protein (target) connection information, and protein (target) information, respectively;
- `path`: str type, the path to store the results, default is `results/`。 If the path does not exist, the corresponding directory will be automatically created.

### `from_TCM_or_Formula`
Start from TCM or formula to construct a TCM integrative pharmacology network.

```python
from DualNet-TCM import analysis
analysis.from_TCM_or_Formula(
  tcm_or_formula, 
  score=990, 
  out_graph=False, 
  out_for_cytoscape=False, 
  out_for_excel=True,
  Research_status_test=False,
  safety_research=False,
  re=True, 
  path='results/'
)
```

`from_TCM_or_Formula` requires a mandatory parameter tcm_or_formula: any combinable data type that can use in to determine if an element is in it, storing the IDs of the TCM or formulas to be analyzed, such as ['DNH0367', 'DNH1695']. Specific ID numbers and TCM/formula information can be found in the [Formula Dataset](/DualNet-TCM/Data/Formula.xlsx)
 and [TCM Dataset](/DualNet-TCM/Data/TCM.xlsx)

`from_TCM_or_Formula`Optional parameters：
- `score`: int type, in the [Chemical_Protein_Links Dataset](/DualNet-TCM/Data/Chemical_Protein_Links.xlsx)
, only records with a combined_score greater than or equal to scorewill be filtered out, default is ；
- `out_graph`: boolean type, whether to output an ECharts-based interactive network visualization graph in HTML format, default is `False`；
- `out_for_cytoscape`: boolean type, whether to output files for Cytoscape visualization, default is `False`；
- `out_for_excel`: boolean  type, whether to output results to an Excel file, default is `True`；
- `research_status_test`: boolean type, whether to evaluate the research value of the obtained targets, default is `False`；
- `safety_research`: boolean type, whether to conduct safety research on the obtained targets, default is `False`；
- `re`: boolean type, whether to return the original analysis results (syndrome differentiation, formula, TCM, compounds (TCM components), proteins (targets), and their connection information), default is True. If re is True, the function will return the results sd, sd_formula_links.xlsx, formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links and proteins, all of which are pd.DataFrame types, storing syndrome differentiation information, syndrome differentiation-formula connection information, formula information, formula-TCM connection information, TCM information, TCM-compound (TCM component) connection information, compound (TCM component) information, compound (TCM component)-protein (target) connection information, and protein (target) information, respectively;
- `path`: str type, the path to store the results, default is `results/`. If the path does not exist, the corresponding directory will be automatically created.

### `from_Protein`
Start from targets to construct TCM integrative pharmacology.

```python
from DualNet-TCM import analysis
analysis.from_Protein(
  protein,
  score=0,
  random_state=None,
  num=1000, 
  tcm_component=True, 
  formula_component=True,
  out_graph=False, 
  out_for_cytoscape=False, 
  out_for_excel=True,
  Research_status_test=False,
  safety_research=False,
  re=True, 
  path='results/'
)
```

`from_Protein` requires a mandatory parameter protein： any combinable data type that can use in to determine if an element is in it, storing the IDs of the TCM or formulas to be analyzed, such as ['ENSP00000381588', 'ENSP00000252519']. Specific ID numbers and target information can be found in the [Target Dataset](/DualNet-TCM/Data/Protein.xlsx)

`from_Protein` Optional parameters:
- `score`: int type, in the [Chemical_Protein_Links Dataset](/DualNet-TCM/Data/Chemical_Protein_Links.xlsx), 
only records with a combined_score greater than or equal to scorewill be filtered out, default is `0`；
- `random_state`: int type, specifies the random seed used by the optimization model, default is `None`, meaning no random seed is specified;
- `num`: int type, specifies the number of solution sets to be generated during optimization, default is `1000`；
- `tcm_component`: boolean type, whether to perform TCM combination optimization, default is `True`；
- `formula_component`: boolean type, whether to perform formula combination optimization, default is `True`；
- `out_graph`: boolean type, whether to output an ECharts-based interactive network visualization graph in HTML format, default is `False`；
- `out_for_cytoscape`: boolean type, whether to output files for Cytoscape visualization, default is `False`；
- `out_for_excel`: boolean type, whether to output results to an Excel file, default is `True`；
- `research_status_test`: boolean type, whether to evaluate the research value of the obtained targets, default is `False`；
- `safety_research`: boolean type, whether to conduct safety research on the obtained targets, default is `False`；
- `re`: boolean type, whether to return the original analysis results (syndrome differentiation, formula, TCM, compounds (TCM components), proteins (targets), and their connection information), default is True. If re is True, the function will return the results sd, sd_formula_links.xlsx, formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links and proteins, all of which arepd.DataFrame  types, storing syndrome differentiation information, syndrome differentiation-formula connection information, formula information, formula-TCM connection information, TCM information, TCM-compound (TCM component) connection information, compound (TCM component) information, compound (TCM component)-protein (target) connection information, and protein (target) information, respectively;
- `path`: str type, the path to store the results, default is `results/`. If the path does not exist, the corresponding directory will be automatically created.

## Results Display

### out_graph

### out_for_cytoscape

### out_for_excel

### research_status_test

### safety_research

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=Carrie-HuYY/DualNet-TCM&type=Date)](https://star-history.com/#Carrie-HuYY/DualNet-TCM&Date)

## Notes

### 1. Regarding Python Version

Since the type annotation syntax like tuple[...] and list[...] is not supported in Python versions before 3.9. Python 3.9 introduced native support for type annotations (PEP 585), but in earlier versions, you need to use the Tuple, List, and other types from the typing module. Therefore, Python ≥3.9 is required. If using a version <3.9, you can modify the compute.py function as follows:

```python
from typing import Tuple, Union

def score(weights: Union[dict, None] = None) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Function logic
    pass
```

### 2.Regarding Data Download
The total size is 9GB. Due to Baidu Netdisk limitations, it is split into three compressed files. After decompression, place them in the data/ folder.

[Download Link 1](https://pan.baidu.com/s/1zIlTjstJMscKdZnP30wc1g?pwd=2n2t) 

[Download Link 2](https://pan.baidu.com/s/1tg8WQtJiJi70A8HIRYG_PA?pwd=9bvh) 

[Download Link 3](https://pan.baidu.com/s/1tg8WQtJiJi70A8HIRYG_PA?pwd=9bvh)
