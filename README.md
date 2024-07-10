# Graph of Frags - GoF tool


Welcome to the `chemicalgof` repository! Graph of Frag (GoF) tool is designed to provide the graph reduction process to convert molecules atom-based to the fragment-based one !
Our tool allows to set a custom-user rule for fragment molecules. By default, our fragmentation rule lead to the so called fragSMILES notation.

To do this, you need python interpreter ... and of course your molecules :)

## Installation
`chemicalgof` package reequires few python dependencies. It currently supports Python 3.11. You can find requirements.txt for dependencies packages.

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
## then get sequence of tokens for fragments and bonds
print(*T.getSequence())

## or simply each fragment and its bonds splitted by dots
print(*T.getString())
```