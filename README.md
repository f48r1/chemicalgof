![python version](https://img.shields.io/badge/python-3.10_|_3.11-white)
![license](https://img.shields.io/badge/license-MIT-orange)
[![Static Badge](https://img.shields.io/badge/ChemRxiv-10.26434/chemrxiv--2024-tm7n6)](https://doi.org/10.26434/chemrxiv-2024-tm7n6)
[![Static Badge](https://img.shields.io/badge/Data%20Zenodo-_10.5281/12713950-blue)](https://doi.org/10.5281/zenodo.12713950)

# Graph of Frags - GoF tool repository

***Molecular Graph Reduction algorithm for fragSMILES notation***

[Introduction](#introduction)
[Installation](#installation)
[How to use](#how-to-use)
[Reference](#reference)

---

## Introduction

Welcome to the `chemicalgof` repository! Graph of Frag (GoF) tool is designed to provide the graph reduction process to convert molecules atom-based to the fragment-based one !
Our tool allows to set a custom-user rule for fragment molecules. By default, our fragmentation rule lead to the so called fragSMILES notation.

To do this, you need python interpreter ... and of course your molecules :)

---

## Installation

`chemicalgof` package reequires few python dependencies. It has been tested through Python>=3.10. You can find requirements.txt for dependencies packages.

The following instructions can be followed to install `chemicalgof` package from our repository:

1. Clone the repository in a comfortable path, for instance `~/`:

    ```shell
    git clone https://github.com/f48r1/chemicalgof.git
    ```

2. (optional but suggested) create your own environment. If you are on a linux OS, python traditional software can be employed. Anaconda is suggested for oher OS like Windows or Mac. For sake of clarity, `gof` is just a customizable name in the following pieces of codes for your virtual environment, therefore you can choose a name by yourself.
    - **python user** traditional python user. Create your own virtual environment for python (and then activate it:

        ```shell
        python -m venv gof
        ```

        ```shell
        source .gof/bin/activate
        ```

    - **conda user** create your own environment for conda (and then activate it):

        ```shell
        conda create --name gof python=3.11
        ```

        `gof` is just a customizable name for your environment, therefore you can choose a name by yourself

        ```shell
        conda activate gof
        ```

3. Navigate to the cloned directory:

    ```shell
    cd chemicalgof/
    ```

4. Install the package :

    ```shell
    python install -m setup.py
    ```

Enjoy :)

---

## How to use

```python
from chemicalgof import Smiles2GoF

## Example SMILES string of a molecule
smiles = 'C[C@@](O)(Cl)C(=O)NC[C@@H]1CC[C@H](C(=O)O)O1' ## molecule provides chirality information

## Convert SMILES to directed graph !
DiG = Smiles2GoF(smiles)
```

Now we need a tokenizer traversing graph to tokenize nodes and edges

```python
from chemicalgof import GoF2Tokens
T=GoF2Tokens(DiG)
```

if you prefer, canonizalied graph can be obtained

```python
from chemicalgof import CanonicalGoF2Tokens
T=CanonicalGoF2Tokens(DiG)
```

```python
### then get sequence of tokens for fragments and bonds
print(*T.getSequence())

## or simply each fragment and its bonds splitted by dots
print(*T.getString())
```

---

## Reference

GoF tool was first presented on the preprint paper [here](https://doi.org/10.26434/chemrxiv-2024-tm7n6)
The resulting fragSMILES notation was used for de novo drug design approaches and compared with traditional notations such as SMILES and SELFIES.

If you think that GoF can be usefull for your project, please cite us :)

```bibtex
@article{Mastrolorito_2024, 
title={fragSMILES: a Chemical String Notation for Advanced Fragment and Chirality Representation}, 
url={http://dx.doi.org/10.26434/chemrxiv-2024-tm7n6}, 
DOI={10.26434/chemrxiv-2024-tm7n6}, 
publisher={American Chemical Society (ACS)}, 
author={Mastrolorito, Fabrizio and Ciriaco, Fulvio and Togo, Maria Vittoria and Gambacorta, Nicola and Trisciuzzi, Daniela and Altomare, Cosimo Damiano and Amoroso, Nicola and Grisoni, Francesca and Nicolotti, Orazio}, 
year={2024}, 
month=jul 
}
```
