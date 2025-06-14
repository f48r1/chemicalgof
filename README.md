![python version](https://img.shields.io/badge/python-3.10_|_3.11-white)
![license](https://img.shields.io/badge/license-MIT-orange)
[![Static Badge](https://img.shields.io/badge/ChemRxiv-10.26434/chemrxiv--2024-tm7n6)](https://doi.org/10.26434/chemrxiv-2024-tm7n6)
[![Static Badge](https://img.shields.io/badge/Data%20Zenodo-_10.5281/12700298-blue)](https://doi.org/10.5281/zenodo.12700298)

> **_NOTE:_**  This package has been refactored and the current version is 0.2.0; Many functions and command series were changed. Deprecated warnings are not implemented yet.

# Graph of Frags - GoF tool repository

**_Molecular Graph Reduction algorithm for fragSMILES notation_**

[Introduction](#introduction)\
[Installation](#installation)\
[How to use](#how-to-use)\
[Reference](#reference)

---

## Introduction

Welcome to the `chemicalgof` repository! Graph of Frag (GoF) tool is designed to provide the graph reduction process to convert molecules atom-based to the fragment-based one !
Our tool allows to set a custom-user rule for fragment molecules. By default, our fragmentation rule lead to the so called fragSMILES notation.

To do this, you need python interpreter ... and of course your molecules :)

---

## 🔧 Installation

### 1. (Optional but recommended) Create a virtual environment

Using a virtual environment is good practice to isolate dependencies.  
You can use either standard Python tools or Conda, depending on your operating system.

- **For Linux users**: the native Python environment is usually sufficient.  
- **For Windows and macOS users**: we recommend using [Anaconda](https://www.anaconda.com/) for better compatibility.


#### 🔹 Using Python `venv`

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

#### 🔹 Using Conda

> ⚠️ In the examples below, `gof` is just a placeholder name for your environment—you can choose any name.

```bash
conda create --name gof python=3.11
conda activate gof
```

---

### 2a. 🔨 Install from source using `setup.py`

1. Clone the repository to a desired directory (e.g., your home folder):

   ```bash
   git clone https://github.com/f48r1/chemicalgof.git
   ```

2. Navigate to the project directory:

   ```bash
   cd chemicalgof/
   ```

3. Install the package locally:

   ```bash
   python setup.py install
   ```

---

### 2b. 📦 Install directly via `pip`

If you prefer a simpler installation, you can install the package directly from GitHub:

```bash
pip install git+https://github.com/f48r1/chemicalgof.git
```


## How to use

```python
from chemicalgof import encode

## Example SMILES string of a molecule
smiles = 'C[C@@](O)(Cl)C(=O)NC[C@@H]1CC[C@H](C(=O)O)O1' ## molecule provides chirality information

## Convert SMILES into relative fragSMILES !
fragsmiles = encode(smiles)

print(fragsmiles)
```

Then, to parse a fragSMILES representation, if valid, and convert it into a molecule

```python
from chemicalgof import decode

ret_smiles = decode(fragsmiles)

# and finally SMILES
print(ret_smiles)
```

Additional detailed examples on how to encode smiles/molecules into fragsmiles are available in [notebook folder](./notebooks/).

---

## Reference

GoF tool was first presented on the preprint paper [here](https://doi.org/10.26434/chemrxiv-2024-tm7n6) and on the actual paper [here](https://www.nature.com/articles/s42004-025-01423-3).
The resulting fragSMILES notation was used for de novo drug design approaches and compared with traditional notations such as SMILES, SELFIES and t-SMILES.

If you think that GoF can be usefull for your project, please cite us :)

```bibtex
@article{mastrolorito_fragsmiles_2025,
title = {{fragSMILES} as a chemical string notation for advanced fragment and chirality representation},
volume = {8},
issn = {2399-3669},
url = {https://doi.org/10.1038/s42004-025-01423-3},
doi = {10.1038/s42004-025-01423-3},
number = {1},
journal = {Communications Chemistry},
author = {Mastrolorito, Fabrizio and Ciriaco, Fulvio and Togo, Maria Vittoria and Gambacorta, Nicola and Trisciuzzi, Daniela and Altomare, Cosimo Damiano and Amoroso, Nicola and Grisoni, Francesca and Nicolotti, Orazio},
month = jan,
year = {2025},
pages = {26},
}
```
