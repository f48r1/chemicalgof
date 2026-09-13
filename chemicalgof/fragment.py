import dataclasses
import weakref
from functools import wraps
from typing import Dict, Tuple

def cached_property(func):
    key = func.__name__

    @property
    @wraps(func)
    def wrapper(self):
        if key not in self._cached:
            self._cached[key] = func(self)
        return self._cached[key]

    return wrapper

from rdkit import Chem

@dataclasses.dataclass(frozen=True, slots=True, weakref_slot=True)
class Fragment:
    smiles: str
    chirality: dict[int, str] = dataclasses.field(default_factory=dict)

    _cached: dict = dataclasses.field(
        init=False,
        repr=False,
        compare=False,
        hash=False,
        default_factory=dict,
    )

    _registry = weakref.WeakValueDictionary()

    def __new__(cls, smiles, chirality=None):
        if chirality is None:
            chirality = dict()
        key = (smiles, tuple(sorted(chirality.items())))

        if (obj := cls._registry.get(key)) is not None:
            return obj

        obj = object.__new__(cls)
        cls._registry[key] = obj
        return obj

    def __getnewargs_ex__(self):
        return (self.smiles,), {"chirality": self.chirality}

    @property
    def mol(self) -> Chem.Mol :
        return Chem.MolFromSmiles(self.smiles)

    @cached_property
    def connector_idxs(self) -> dict[int, int] :
        return {
            a.GetIdx() : num_Hs
            for a in self.mol.GetAtoms()
            if (num_Hs := a.GetTotalNumHs()) > 0
        }

    @property
    def num_connector(self) -> int:
        return len(self.connector_idxs)

    @cached_property
    def group_symmetric_connector(self) -> tuple[tuple[int]] :

        symm_dict = {}

        for atom_idx, order in enumerate(Chem.CanonicalRankAtoms(self.mol, breakTies=False, includeChirality=False)):

            if atom_idx not in self.connector_idxs:
                continue
            elif order not in symm_dict:
                symm_dict[order] = []

            symm_dict[order].append(atom_idx)

        return tuple([ tuple(group) for group in symm_dict.values() ])
    
    @property
    def suffix(self) -> str:
        if not self.chirality:
            return ''
        elif self.num_connector == 1:
            return list(self.chirality.values())[0] # Only one chiral atom does not include its index (e.g. C|S)
        else:
            return ''.join([str(k)+v for k,v in self.chirality.items()])

    @property
    def fragsmiles(self) -> str:
        return self.smiles + ('|' + self.suffix if self.suffix else '')
    
    def __repr__(self):
        return self.fragsmiles
    
    @cached_property
    def typeId(self):

        # FIXME by including stereo or not ?
        
        # NOTE string representation by not including stereochemistry ...
        string = self.smiles

        # NOTE string representation includes stereochemistry ;)
        # string = self.fragsmiles
        
        return int("".join([str(ord(c)) for c in string]))
