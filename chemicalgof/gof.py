from .utils import GetPotAtomLinkers
from typing import Union

import networkx as nx, numpy as np
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Classes for reduced graph handling

def fNodeMatch(n1,n2):
    return all( [
        n1["smiles"]==n2["smiles"],
        n1.get("chirality")==n2.get("chirality"), # TODO check me
    ])

def fEdgeMatch(e1,e2):
    return all([
        e1["aB"]==e2["aB"],
        e1.get('chirality')==e2.get('chirality'), # FIXME
               ])

# related class for single frag in reduced graph
class FragNode:

    _IDs = { } # NOTE this private attributes need just to store every unique integer number relative to different smiles... it is not very required but useful.

    def __init__(self, smiles:str|None = None, chirality: dict[int,str] | None = None, parent : Union['DiGraphFrags', None] = None):
        # NOTE parent is set when node is added to Graph at the first time.
        self.parent = parent

        if chirality is not None:
            self.chirality = chirality
        else:
            self.chirality = {}

        if smiles is not None:
            self.smiles = smiles
            self.PotAtomLinkers:list[int]=GetPotAtomLinkers(smiles)
        else:
            self.smiles = ''
            self.PotAtomLinkers=[]

    def __repr__(self):
        if self.idx>=0:
            return f"{self.smiles} idx:{self.idx}"
        else:
            return f"{self.smiles}"
        
    @property
    def suffix(self) -> str:
        if not self.chirality:
            return ''
        elif self.numPotAtomLinkers == 1:
            return list(self.chirality.values())[0] # Only one chiral atom
        else:
            return ''.join([str(k)+v for k,v in sorted(self.chirality.items())])
        
    def __str__(self):
        return self.smiles + ('|' + self.suffix if self.suffix else '')

    def setParent(self, parent):
        self.parent=parent
    
    @property
    def idx(self):
        if self.parent:
            return list(self.parent.nodes).index(self)
        else:
            return -1
        
    @property
    def typeId(self):
        
        id:int = self._IDs.get(self.smiles)
        if id is None:
            # id = self._IDs[self.smiles] = int("".join([str(ord(c)) for c in self.__str__() ])) # NOTE string representation includes stereochemistry ;)
            id = self._IDs[self.smiles] = int("".join([str(ord(c)) for c in self.smiles])) # NOTE string representation by not including stereochemistry ...
        return id
    
    @property
    def numPotAtomLinkers(self):
        return len(self.PotAtomLinkers)

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
        
    def _set_details(self, node: FragNode):
        if not node.parent:
            node.setParent(self)
        if "smiles" not in nx.get_node_attributes(self, node):
            self.nodes[node]["smiles"]=node.smiles
        if "chirality" not in nx.get_node_attributes(self, node):
            self.nodes[node]["chirality"]=node.chirality

    def add_node(self, node_for_adding: FragNode, **attr):
        nx.DiGraph.add_node(self, node_for_adding, **attr)
        self._set_details(node_for_adding)
        
    def add_nodes_from(self, nodes_for_adding, **attr):
        nx.DiGraph.add_nodes_from(self,nodes_for_adding, **attr )
        for node in nodes_for_adding:
            if np.iterable(node):
                if type(node) != dict:
                    node=node[0]
            self._set_details(node)
