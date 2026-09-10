from rdkit import Chem
import warnings

LEGACY_LABEL_IS_CORRECT = True

# [ ] something to do here ?

def GetPotAtomLinkers(s):
    if "|" in s: # if fragment has chirality
        s, *_=s.split("|")
    m=Chem.MolFromSmiles(s)
    
    return [a.GetIdx() for a in m.GetAtoms() if a.GetTotalNumHs()>0]

def CanonizeMol(m, retMap=False):
    for a in m.GetAtoms():
        a.SetProp("oldIdx", str(a.GetIdx()))
    _ = Chem.MolToSmiles(m)
    
    order = m.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"]
    
    m_canonical = Chem.RenumberAtoms(m, order)
    mapIdx={int(a.GetProp("oldIdx")) : a.GetIdx() for a in m_canonical.GetAtoms()}
    
    if retMap:
        return m_canonical, mapIdx
    else:
        return m_canonical

def ClearSmiles(s):
    m=Chem.MolFromSmiles(s)
    m=Chem.ReplaceSubstructs(m, Chem.MolFromSmiles("*"), Chem.MolFromSmiles("[H]"), replaceAll=True)[0]
    # m=Chem.AddHs(m)
    m=Chem.RemoveAllHs(m)
    
    return Chem.MolToSmiles(m)

def CanonizeFragWithDummies(m1):
    # m1=Chem.MolFromSmiles(smi)
    m2,aIdxsMap1=CanonizeMol(m1,True)
    aIdxs2=[a.GetIdx() for a in m2.GetAtoms() if a.GetAtomicNum()!=0 ]
    m3=Chem.RemoveAllHs(  Chem.ReplaceSubstructs(m2, Chem.MolFromSmiles("*"), Chem.MolFromSmiles("[H]"), replaceAll=True)[0]  )
    aIdxs3=[a.GetIdx() for a in m3.GetAtoms() if a.GetAtomicNum()!=0 ]
    
    tmp={ k:aIdxs3[aIdxs2.index(v)] for k,v in aIdxsMap1.items() if v in aIdxs2}
    if Chem.MolToSmiles(m3)!=Chem.MolToSmiles(m3, canonical=False):
        m4, aIdxsMap3=CanonizeMol(m3, True)
        retMap={ k:aIdxsMap3[v] for k,v in tmp.items() if v in aIdxs3}
        return m4, retMap
    
    return m3, tmp


def FindProperStereoCenters(mol:Chem.Mol, warning=True) -> dict[int, str]:

    with_legacy : dict[int,str] = dict( Chem.FindMolChiralCenters(mol, useLegacyImplementation=True ) )

    # XXX NOTE useLegacyImplementation=False -> r,s CipLabel also included
    try:
        without_legacy : dict[int,str] = dict( Chem.FindMolChiralCenters(mol, useLegacyImplementation=False ) )
    except Exception as e:
        warnings.warn(
            f"Could not assign modern CIP labels: \n{e}",
            RuntimeWarning,
        )
        without_legacy = {}

    # NOTE this order matters: without_legacy will replace psuedo chirality (r or s)
    if without_legacy and with_legacy:

        corrupted_cip_labels = {}

        for atom_idx in tuple(with_legacy.keys()):

            cip_label = with_legacy.pop(atom_idx)
            prev_cip_label = without_legacy.get(atom_idx, None)
            
            if prev_cip_label is None:
                without_legacy[atom_idx] = cip_label
            elif prev_cip_label.upper() != cip_label.upper():
                corrupted_cip_labels[atom_idx] = cip_label

        if not corrupted_cip_labels:
            return without_legacy
        elif not warning:
            if LEGACY_LABEL_IS_CORRECT:
                return without_legacy | corrupted_cip_labels
            else:
                return without_legacy

        warning_msg = '\n' + 'Opposite CIP labels from different Legacy implementations.' + '\n' + 'CIP labels '
        if LEGACY_LABEL_IS_CORRECT:
            warning_msg += 'from Legacy implementation '
            labels_str = ', '.join([f'{atom_idx}{cip_label}' for atom_idx, cip_label in corrupted_cip_labels.items() ])
            without_legacy.update(corrupted_cip_labels)
        else:
            warning_msg += 'without Legacy implementation '
            labels_str = ', '.join([f'{atom_idx}{without_legacy[atom_idx]}' for atom_idx in corrupted_cip_labels.keys() ])

        corrupted_cip_labels.clear()

        warning_msg += 'were employed:' + '\n' + labels_str
        warnings.warn(warning_msg)

        return without_legacy
    
    elif not without_legacy and not with_legacy:
        return {}
    elif without_legacy and not with_legacy:
        return without_legacy
    elif not without_legacy and with_legacy:
        return with_legacy
