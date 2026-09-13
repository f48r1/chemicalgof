import networkx as nx

# Classe for nodes
from .fragment import Fragment
import dataclasses

@dataclasses.dataclass(frozen=True, eq=False, slots=True)
# related class for single frag in reduced graph
class FragNode:
    fragment: Fragment

def fEdgeMatch(e1,e2):
    return all([
        e1["aB"]==e2["aB"],
        e1.get('stereo')==e2.get('stereo'),
               ])

# related class for directed reduced graph
class DiGraphFrags(nx.DiGraph):

    def __init__(self):
        super().__init__()
        
    def __eq__(self, other):
        if type(other) is DiGraphFrags:
            # TODO does equality works in this way ?
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
