from rdkit import Chem
import numpy as np
import pandas as pd
from typing import Union, Literal

from .exceptions import InvalidChirality
from .gof import DiGraphFrags

# return RDKit mol. If invalid reduced graph, it raises error
class Assembler:

    attempts_multiplier = 10 # NOTE rly important ! Maybe ...

    def __init__(self, DiG:DiGraphFrags, strict_chirality:bool = True):
        self.DiG = DiG
        self.stereo_atoms = {}
        self.strict_chirality = strict_chirality

    def assemble_mol(self):
        # Setup indexes
        frag_smiles=[_.smiles for _ in self.DiG._node]
        frag_num_atoms=[Chem.MolFromSmiles(_).GetNumAtoms() for _ in frag_smiles]

        additive_idx=frag_num_atoms[:-1]
        additive_idx=np.cumsum(additive_idx).tolist()
        additive_idx.insert(0,0)

        additive_idx={k:v for k,v in zip(self.DiG._node, additive_idx)}
        frag_num_atoms={k:v for k,v in zip(self.DiG._node,frag_num_atoms)}

        # Pieces of molecule to be edited
        editable_mol=Chem.EditableMol( Chem.MolFromSmiles(".".join(frag_smiles)) )
        disconnected_mol=editable_mol.GetMol() # Raw molecule without bonds. need for Hydrogens count

        # Chirality from node attributes
        self.stereo_atoms.update(
            {
                atom+additive_idx[n]:stereo
                for n in self.DiG._node if n.chirality
                for atom,stereo in n.chirality.items()
            }
        )
        
        # Setting edges
        edges=list(self.DiG.edges)
        explHs={}
        while edges:
            a,b = edges[0]
            datas= ( self.DiG.get_edge_data(a, b), self.DiG.get_edge_data(b,a) )
            if datas[0]["aB"] >= frag_num_atoms[a] or datas[1]["aB"] >= frag_num_atoms[b]:
                raise ValueError('Bond Error : Invalid atom index for a fragment')
            bond = datas[0]["aB"]+additive_idx[a], datas[1]["aB"]+additive_idx[b]
            editable_mol.AddBond(*bond, order=Chem.rdchem.BondType.SINGLE)
            edges.remove((a,b))
            edges.remove((b,a))

            # Setting explicited Hydrogens and chirality from edges
            for data, atom in zip(datas, bond ):
                if data['stereo'] is not None:
                    self.stereo_atoms[atom]=data["stereo"]
                if disconnected_mol.GetAtomWithIdx(atom).GetNumExplicitHs():
                    if atom in explHs:
                        explHs[atom]+=1
                    else:
                        explHs[atom]=1

        gen=editable_mol.GetMol()

        gen = Chem.RemoveAllHs(gen, sanitize=False)
        for a,Hs in explHs.items():
            atom=gen.GetAtomWithIdx(a)
            numExplHs = atom.GetNumExplicitHs()
            if Hs > numExplHs:
                raise ValueError(f'Explicit valence for atom # {a} {atom.GetSymbol()} is greater than permitted')
            atom.SetNumExplicitHs(numExplHs-Hs)
        Chem.SanitizeMol(gen)

        self.mol:Chem.Mol = gen
    
    def initialize_chirality(self):
        self.chir_table = chir_table = pd.DataFrame.from_dict(self.stereo_atoms, columns=['tgt'], orient='index')
        chir_table['current'] = ''

        # Initializing stereochemistry tag
        for atom_idx in self.stereo_atoms.keys():
            chirA=self.mol.GetAtomWithIdx(atom_idx)
            chirA.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

    def update_chirality(self, check:bool=True, cleanIt:bool=False):
        Chem.AssignStereochemistry(self.mol, cleanIt=cleanIt, force=True)
        Chem.AssignCIPLabels(self.mol, list(self.stereo_atoms.keys()))  

        for atom_idx, cip_label in self.stereo_atoms.items():
            chirAtom=self.mol.GetAtomWithIdx(atom_idx)
            if chirAtom.HasProp("_CIPCode") == 0 and self.strict_chirality and check:
                raise InvalidChirality(chirAtom.GetSymbol(), atom_idx, cip_label)
            elif chirAtom.HasProp("_CIPCode") == 0:
                self.chir_table.at[atom_idx,'current'] = ''
            else:
                self.chir_table.at[atom_idx,'current'] = chirAtom.GetProp("_CIPCode")

    def get_mask_wrong(self, _type=['R','S']):
        return (
            (self.chir_table['tgt'].str.upper() != self.chir_table['current'].str.upper()) &
            (self.chir_table['current'].notnull()) &
            (self.chir_table['current'].isin(_type))
        )

    def invert_wrong_chirality(self, _type: Union[Literal['lower'] ,Literal['capital']] = 'capital'):
        typ=['r','s'] if _type == 'lower' else ['R','S']

        mask_wrong = self.get_mask_wrong(typ)
        wrongs = mask_wrong.sum()
        
        attempts_max = self.chir_table['current'].isin(typ).sum() * self.attempts_multiplier
        attempts = 0
        while wrongs > 0 and attempts < attempts_max:
            row = self.chir_table[mask_wrong].sample(n=1).squeeze()

            atom_idx = int(row.name)
            atom = self.mol.GetAtomWithIdx(atom_idx)
            atom.InvertChirality()
            if _type == 'capital':
                self.update_chirality(check=False, cleanIt=False)
            else:
                self.update_chirality()

            mask_wrong = self.get_mask_wrong(typ)
            wrongs = mask_wrong.sum()
            attempts+=1

            if _type != 'lower':
                continue

            mask_capital_wrong = self.get_mask_wrong(['R','S'])
            if mask_capital_wrong.sum () > 0:
                if self.strict_chirality:
                    raise InvalidChirality()
                atom.InvertChirality()
                self.update_chirality()
                atom.SetChiralTag(Chem.CHI_UNSPECIFIED)
                self.update_chirality()

def GoF2Mol(DiG:DiGraphFrags, strict_chirality:bool=True) -> 'Chem.Mol':
    """Explode reduced graph into rdkit mol object.

    Args:
        DiG (DiGraphFrags): Reduced graph
        strict_chirality (bool, optional): If consider invalid chirality labels provided for atoms. Defaults to True.

    Raises:
        InvalidChirality: If scrict_chirality is True and invalid chiality labels are provided.

    Returns:
        Chem.Mol: RDKit mol object obtained by exploding reduced graph.
    """

    assembler = Assembler(DiG, strict_chirality)
    assembler.assemble_mol()

    if not assembler.stereo_atoms:
        return assembler.mol

    assembler.initialize_chirality()
    assembler.update_chirality(check=False)

    assembler.invert_wrong_chirality('capital')

    assembler.invert_wrong_chirality('lower')

    assembler.update_chirality( check=True, cleanIt=True)

    if strict_chirality:
        current_stereo = dict(Chem.FindMolChiralCenters(assembler.mol, useLegacyImplementation=False)) # BUG we need either legacy=True and False
        current_stereo = current_stereo | dict(Chem.FindMolChiralCenters(assembler.mol, useLegacyImplementation=True))

        for atom_idx, cip_label in assembler.stereo_atoms.items():
            current_label = current_stereo.get(atom_idx)
            chirAtom = assembler.mol.GetAtomWithIdx(atom_idx)
            if current_label is None or current_label.upper() != cip_label.upper():
                raise InvalidChirality(chirAtom.GetSymbol(), atom_idx, cip_label)

    return assembler.mol
