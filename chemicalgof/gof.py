from typing import Union

import networkx as nx, numpy as np
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Classes for reduced graph handling


from .fragment import Fragment
import dataclasses

@dataclasses.dataclass(frozen=True, eq=False, slots=True)
# related class for single frag in reduced graph
class FragNode:
    fragment: Fragment

def fEdgeMatch(e1,e2):
    return all([
        e1["aB"]==e2["aB"],
        e1.get('stereo')==e2.get('stereo'), # FIXME
               ])

# related class for directed reduced graph
class DiGraphFrags(nx.DiGraph):

    def __init__(self):
        super().__init__()
        
    def __eq__(self, other):
        if type(other) is DiGraphFrags:
            return nx.is_isomorphic(self, other, edge_match=fEdgeMatch)
        return False
    
    def GetFragsByIdx(self,*idxs:int) -> list[int] | int:
        nodes = [node for idx, node in enumerate(self._node.keys()) if idx in idxs]
        if len(idxs) == 1:
            return nodes[0]
        return nodes

    def GetEdgesByIdx(self,*idxs):
        edges = [edge for idx, edge in enumerate(self.edges.keys()) if idx in idxs]
        if len(idxs) == 1:
            return edges[0]
        return edges
        
    # def _set_details(self, node: FragNode):
    #     if not node.parent:
    #         node.setParent(self)
    #     if "smiles" not in nx.get_node_attributes(self, node):
    #         self.nodes[node]["smiles"]=node.smiles
    #     if "chirality" not in nx.get_node_attributes(self, node):
    #         self.nodes[node]["chirality"]=node.chirality

    # def add_node(self, node_for_adding: FragNode, **attr):
    #     nx.DiGraph.add_node(self, node_for_adding, **attr)
    #     self._set_details(node_for_adding)
        
    # def add_nodes_from(self, nodes_for_adding, **attr):
    #     nx.DiGraph.add_nodes_from(self,nodes_for_adding, **attr )
    #     for node in nodes_for_adding:
    #         if np.iterable(node):
    #             if type(node) != dict:
    #                 node=node[0]
    #         self._set_details(node)
