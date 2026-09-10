from rdkit.Chem import Draw
from rdkit import Chem

from ..gof import FragNode

def set_options(d2d):
    dopts = d2d.drawOptions()

    dopts.addAtomIndices = True
    dopts.noAtomLabels = False
    dopts.explicitMethyl = True
    dopts.annotationFontScale = 1.0

def drawMol(mol:Chem.RWMol, as_svg=False):

    if as_svg:
        d2d = Draw.MolDraw2DSVG(-1,-1)
    else:
        d2d = Draw.MolDraw2DCairo(-1,-1)

    set_options(d2d)

    for atom in mol.GetAtoms():
        atom:Chem.Atom
        # atom.SetProp("atomLabel", atom.GetSymbol())
        atom.SetNoImplicit(True)
    # mol.UpdatePropertyCache()

    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()

    data_img = d2d.GetDrawingText()

    if not as_svg:
        from io import BytesIO
        from PIL import Image
        bio = BytesIO(data_img)
        return Image.open(bio)

    data_img = data_img.replace('svg:','')
    
    import re
    data_img = re.sub(r'(?m)^<rect.*\n', '', data_img)

    return data_img

def drawNode(node:FragNode, as_svg=False):
    if as_svg:
        d2d = Draw.MolDraw2DSVG(-1,-1)
    else:
        d2d = Draw.MolDraw2DCairo(-1,-1)

    dopts = d2d.drawOptions()
    dopts.noAtomLabels = False

    fragment = node.fragment
    mol = Chem.MolFromSmiles(fragment.smiles)

    def idxs2atoms(atom_idxs):
        return [mol.GetAtomWithIdx(atom_idx) for atom_idx in atom_idxs]

    def CommonAtoms(atom_idxs):
        for atom in idxs2atoms(atom_idxs):
            atom.SetNoImplicit(True)

    def Terminal_sp2_Cs(atom_idxs):
        for atom in idxs2atoms(atom_idxs):
            atom.SetProp("atomLabel", 'C')
            atom.SetNoImplicit(True)

    def ChiralAtoms(chirality_dict:dict[int,str]):
        for atom_idx, cip_label in chirality_dict.items():
            atom = mol.GetAtomWithIdx(atom_idx)
            atom.SetProp('atomNote', cip_label)
            atom.SetNoImplicit(True)

    def FixSingleStereoAtom(atom_idx, cip_label):
        single_stereo_atom = mol.GetAtomWithIdx(single_stereo_atom_idx)

        single_stereo_atom.SetNoImplicit(True)
        single_stereo_atom.SetProp("atomLabel", single_stereo_atom.GetSymbol() + '|' + single_stereo_label)

    terminal_sp2_Cs = tuple( match[0] for match in mol.GetSubstructMatches(Chem.MolFromSmarts('[C^2H2]')) )
    all_atom_idxs = set(range(mol.GetNumAtoms()))

    if fragment.num_connector > 1:
        dopts.addAtomIndices = True
        dopts.annotationFontScale = 1.0
        
    elif fragment.num_connector == 1 and len(fragment.chirality) == 1:

        (single_stereo_atom_idx, single_stereo_label),*_ = fragment.chirality.items()
        FixSingleStereoAtom(single_stereo_atom_idx, single_stereo_label)
        all_atom_idxs.discard(single_stereo_atom_idx)

    elif fragment.chirality:
        ChiralAtoms(fragment.chirality)
        all_atom_idxs.difference_update(fragment.chirality.keys())

    if terminal_sp2_Cs:
        Terminal_sp2_Cs(terminal_sp2_Cs)
        all_atom_idxs.difference_update(terminal_sp2_Cs)

    CommonAtoms(all_atom_idxs)

    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()

    data_img = d2d.GetDrawingText()

    if not as_svg:
        from io import BytesIO
        from PIL import Image
        bio = BytesIO(data_img)
        return Image.open(bio)

    data_img = data_img.replace('svg:','')
    
    import re
    data_img = re.sub(r'(?m)^<rect.*\n', '', data_img)

    return data_img