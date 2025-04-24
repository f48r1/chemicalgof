from rdkit import Chem

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


