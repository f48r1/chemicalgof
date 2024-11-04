from .utils import GetPotAtomLinkers
import networkx as nx, numpy as np, re
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Classes for reduced graph handling

def fNodeMatch(n1,n2):
    return all( [
        n1["smiles"]==n2["smiles"],
        n1.get("suffStereo")==n2.get("suffStereo"),
    ])

def fEdgeMatch(e1,e2):
    return all([
        e1["aB"]==e2["aB"],
        e1.get('stereo')==e2.get('stereo'),
               ])

# related class for single frag in reduced graph
class FragNode:
    @classmethod
    def fromSmiles(cls, smiles):
        ret=cls()

        if "|" in smiles:
            # if smiles provides stereochemical information, smiles and info are splitted
            ret.smiles, ret.suffStereo=smiles.split("|")
        else:
            ret.smiles=smiles; ret.suffStereo=None

        # fixed data related to frag
        ret.PotAtomLinkers=GetPotAtomLinkers(ret.smiles)
        ret.numPotAtomLinkers=len(ret.PotAtomLinkers)

        ## typeId is unique integer number relative to the smiles string. We need it to the canonicalization process
        ret.typeId = int("".join([str(ord(c)) for c in smiles]))
        
        return ret

    def __init__(self):
        ### FORSE VOGLIO FISSARE IL PARENT QUANDO VIENE AGGIUNTO PER LA PRIMA VOLTA AD UN NODO. COSI RIMANE IL SUO E DI NESSUN'ALTRO
        self.parent=None        
        
    def __repr__(self):
        if self.idx>=0:
            return f"{self.smiles} idx:{self.idx}"
        else:
            return f"{self.smiles}"

    def setParent(self, parent):
        self.parent=parent
    
    @property
    def idx(self):
        if self.parent:
            return list(self.parent.nodes).index(self)
        else:
            return -1

# related class for undirected reduced graph
class GraphFrags(nx.Graph):
    def __init__(self):
        super().__init__()
 
    def __eq__(self, other):
        if type(other) is GraphFrags:
            return nx.is_isomorphic(self,other, node_match=fNodeMatch)
        return False
    
    def GetFragsByIdx(self, *args):
        self.GetFragsByIdx = self.refG.GetFragsByIdx
        return self.GetFragsByIdx(*args)
    
    def GetEdgesByIdx(self, *args):
        self.GetEdgesByIdx = self.refG.GetEdgesByIdx
        return self.GetEdgesByIdx(*args)

# related class for directed reduced graph
class DiGraphFrags(nx.DiGraph):
    def __init__(self):
        super().__init__()
        
    def __eq__(self, other):
        if type(other) is DiGraphFrags:
            return nx.is_isomorphic(self,other, node_match=fNodeMatch, edge_match=fEdgeMatch)
        return False
    
    def to_undirected_class(self):
        return GraphFrags
    
    def to_undirected(self):
        G=nx.Graph.to_undirected(self)
        G.refG=self
        return G
    
    def GetFragsByIdx(self,*args):
        if len(list(args))>1:
            return np.array(self.nodes)[list(args)]  
        else: 
            return np.array(self.nodes)[list(args)][0]

    def GetEdgesByIdx(self,*args):
        if len(list(args))>1:
            return np.array(self.edges)[list(args)]  
        else: 
            return np.array(self.edges)[list(args)][0]
        
    def add_node(self, node_for_adding, **attr):
        nx.DiGraph.add_node(self, node_for_adding, **attr)
        if not node_for_adding.parent:
            node_for_adding.setParent(self)
        if "smiles" not in nx.get_node_attributes(self, node_for_adding):
            self.nodes[node_for_adding]["smiles"]=node_for_adding.smiles
        if "suffStereo" not in nx.get_node_attributes(self, node_for_adding):
            self.nodes[node_for_adding]["suffStereo"]=node_for_adding.suffStereo
        
    def add_nodes_from(self, nodes_for_adding, **attr):
        nx.DiGraph.add_nodes_from(self,nodes_for_adding, **attr )
        for node in nodes_for_adding:
            if np.iterable(node):
                if type(node) != dict:
                    node=node[0]
            if not node.parent:
                node.setParent(self)
            if "smiles" not in nx.get_node_attributes(self, node):
                self.nodes[node]["smiles"]=node.smiles
            if "suffStereo" not in nx.get_node_attributes(self, node):
                self.nodes[node]["suffStereo"]=node.suffStereo

    # return RDKit mol. If invalid reduced graph, it raises error
    def getMol(self):
        smis=[_.smiles for _ in self._node]
        totAtom=[Chem.MolFromSmiles(_).GetNumAtoms() for _ in smis]
        cumAtom=totAtom[:-1]
        cumAtom=np.cumsum(cumAtom).tolist()
        cumAtom.insert(0,0)
        cumAtom={k:v for k,v in zip(self._node, cumAtom)}
        totAtom={k:v for k,v in zip(self._node,totAtom)}

        editMol=Chem.EditableMol( Chem.MolFromSmiles(".".join(smis)) )
        _tmpMol=editMol.GetMol()

        edges=list(self.edges)
        explHs={}
        stereoAtoms={ int(s[:-1])+cumAtom[n]:s[-1]  for n in self._node if n.suffStereo!=None 
                         for s in re.findall(r"\d+[A-Z]{1}",n.suffStereo) }
        
        while edges:
            a,b = edges[0]
            datas= ( self.get_edge_data(a, b), self.get_edge_data(b,a) )
            if datas[0]["aB"] >= totAtom[a] or datas[1]["aB"] >= totAtom[b]:
                raise Exception('invalid atom index for a fragment')
            bond = datas[0]["aB"]+cumAtom[a], datas[1]["aB"]+cumAtom[b]
            editMol.AddBond(*bond, order=Chem.rdchem.BondType.SINGLE)
            edges.remove((a,b))
            edges.remove((b,a))

            for data, atom in zip(datas, bond ):
                if "stereo" in data: #and not atom in stereoAtoms:
                    stereoAtoms[atom]=data["stereo"]
                if _tmpMol.GetAtomWithIdx(atom).GetNumExplicitHs():
                    if atom in explHs:
                        explHs[atom]+=1
                    else:
                        explHs[atom]=1

        gen=editMol.GetMol()

        # print(explHs)
        # return gen
        gen = Chem.RemoveAllHs(gen, sanitize=False)
        for a,Hs in explHs.items():
            atom=gen.GetAtomWithIdx(a)
            atom.SetNumExplicitHs(atom.GetNumExplicitHs()-Hs)
        Chem.SanitizeMol(gen)

        for k,v in stereoAtoms.items():
            chirA=gen.GetAtomWithIdx(k)
            chirA.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            Chem.AssignStereochemistry(gen, cleanIt=True, force=True)
            if chirA.GetProp("_CIPCode").upper() != v.upper():
                chirA.InvertChirality()
                Chem.AssignStereochemistry(gen, cleanIt=True, force=True)

        return gen