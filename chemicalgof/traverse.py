import itertools
import numpy as np
import networkx as nx
from .gof import FragNode, DiGraphFrags, GraphFrags
import random
secure_random = random.SystemRandom()

# NOTE This is like a pointer class of FragNode which has additional attributes to interact with branching paths
class Pathstep:
    def __init__(self, frag_node:FragNode):
        self.frag_node = frag_node # pointing to the node
        self.branches: list[list[Pathstep]]=[]

    def __repr__(self):
        b=f"({len(self.branches)})" if self.branches else ""
        return str(self.frag_node)+b


class Traverser:
    def __init__(self, DiG: DiGraphFrags, canonize:bool = True, casual: bool | int = False):
        self.DiG = DiG
        self.G = DiG.to_undirected()

        self.canonize = canonize
        self.casual = casual

    @staticmethod
    def getUniqueMinSum(paths: list[list[Pathstep]]) -> int | None:
        """
        Integer ID number identifiying FragNodes are summed sequentially for each path. Pathway index relative to the lowest sum value is returned.
        Args:
            paths (list[list[Pathstep]]): Pathways to be compared by FragNodes ID sum.

        Returns:
            int | None: Pathway index of input list returning the lowest sum value. None if equal values are recognized.
        """
        s=[0 for _ in paths ]
        for steps in zip(*paths):
            for idx, step in enumerate(steps):
                s[idx] += step.frag_node.typeId
            minimum = min(s)
            _min = [pos for pos,element in enumerate(s) if element==minimum]
            if len(_min)==1 :
                return _min[0]
        return None ## No any path ranked as first

    @staticmethod
    def getMultiMinSum(paths: list[list[Pathstep]]) -> list[int]:
        """
        Integer ID number identifiying FragNodes are summed sequentially for each path. Pathway indexes relative to the lowest sum value are returned.
        Args:
            paths (list[list[Pathstep]]): _description_

        Returns:
            list[int]: Pathway indexes of input list returning the lowest sum value.
        """
        s=[0 for _ in paths ]
        idxs=np.array(range(len(paths)))
        for steps in zip(*paths):
            IDs=[step.frag_node.typeId for idx,step in enumerate(steps) if idx in idxs]

            s = [a+b  for a,b in zip (s,IDs)]
            minimum = min(s)
            _min = [pos for pos,element in enumerate(s) if element==minimum]
            if len(_min)==1 :
                return _min
            elif len(_min)<len(idxs):
                idxs = idxs[_min]
                s = [element for i,element in enumerate(s) if i in _min ]
        return idxs

    @staticmethod
    def groupFragsByType(paths : list[list[Pathstep]]) -> list[list[int]]:
        """_summary_

        Args:
            paths (list[list[Pathstep]]): _description_

        Returns:
            list[list[int]]: _description_
        """
        pathsFrags=[ tuple(str(step.frag_node.typeId) for step in path) for path in paths ]
        groups=[ sorted(np.where( (np.array(pathsFrags)==k).all(axis=1))[0] ) 
                for k in set(pathsFrags) ]
        return groups
    
    def find_paths(self, source: FragNode | None = None, subG:GraphFrags | None = None,) -> list[list[Pathstep]] :
        G = self.G if subG is None else subG

        terminalNodes = [node for node in G._node if G.degree[node]==1 and node != source]
        # if the source of the path is not specified, all paths between terminal nodes are recognized
        if not source:
            terminalNodesCouple=[[(a,b),(b,a)] for a,b in itertools.combinations(terminalNodes,2)]
            terminalNodesCouple=itertools.chain(*terminalNodesCouple)
        # only paths from source to terminal nodes are recognized
        else:
            terminalNodesCouple=itertools.product([source],terminalNodes)

        ## path composed of 1 node ... :)
        if not terminalNodes and source is not None:
            return [[Pathstep(source)]]
        elif not terminalNodes and source is None:
            return None # BUG rdkit does not recognize node !

        paths=[ list(nx.all_simple_paths(G, a,b))[0] for a,b in terminalNodesCouple ]

        return [ [Pathstep(step) for step in path] for path in paths ]
    
    def filter_paths(self, paths:list[list[Pathstep]], lenMax:int=-1) -> list[list[Pathstep]] :
        current = 0 if lenMax < 1 else lenMax
        selected = []
        for path in paths:
            l=len(path)
            if lenMax < 1 and l > current: # NOTE if we are looking for longest paths
                selected=[path]
                current=l
            if l==current:
                selected.append(path)

        return selected

    # path input is referred to the main path of reduced graph. So first and last element of list 'path' are terminal nodes
    def findFirstBranchIdx(self, path:list[Pathstep]) -> int | None:
        # NOTE we can ignore terminal nodes cause of they have just 1 deegree for sure.
        for i,step in enumerate(path[1:-1]):
            node = step.frag_node
            if self.G.degree[node]>2:
                return i+1
        return np.nan ## it is linear.
    
    def traverse_by_idx(self, nodes:list[FragNode] | None = None, first_node:FragNode | None = None) -> list[Pathstep]:

        if nodes is None:
            nodes = list(self.DiG._node)
        else:
            nodes.sort(key = lambda x: x.idx)

        if first_node is None:
            for node in nodes:
                if self.G.degree[node]==1:
                    first_node = node
                    break

        pathway = [first_node]
        nodes.remove(first_node)

        i=0
        while i < len(nodes):

            next_node = nodes[i]
            last_node = pathway[-1]

            if next_node in self.G.neighbors(last_node):
                nodes.remove(next_node)
                pathway.append(next_node)
                i=0

                if self.G.degree[next_node]==1:
                    break
            else:
                i+=1
            
        pathway = [Pathstep(node) for node in pathway]

        return pathway

    # Recursive function to build nested pathways
    def buildBranches(self, main_path:list[Pathstep], ascendent:FragNode | None = None) -> list[Pathstep]:
        G = self.G

        main_path_nodes = [step.frag_node for step in main_path]

        for step in main_path:
            frag_node = step.frag_node
            branching_initilizer_nodes: list[FragNode]=[neighbor for neighbor in G.neighbors(frag_node) 
                    if neighbor not in main_path_nodes and neighbor!=ascendent]
            
            if not branching_initilizer_nodes:
                continue

            branching_paths = []
            for branching_node in branching_initilizer_nodes:
                branching_graph=self.split_isolate_between(frag_node, branching_node)
                branching_pathway=self.findBranchPath(branching_graph, branching_node)
                self.buildBranches(branching_pathway, ascendent=frag_node)

                branching_paths.append(branching_pathway)

            idxSorted = self.sortBranches(branching_paths)
            for i in idxSorted:
                step.branches.append(branching_paths[i])

        return main_path
    
    def split_isolate_between(self, core:FragNode, source:FragNode) -> GraphFrags:
        # NOTE core is the initializer node of branching path
        subG=self.G.copy()
        subG.remove_edge(core, source)
        toDel=[x for x in nx.connected_components(subG) if source not in x]
        subG.remove_nodes_from(*toDel)
        return subG # [x] we could return a subgraph but it's the same ...

    def getNumBranches(self, path, subG=None):
        G = self.G if subG is None else subG
        return sum([G.degree[step.frag_node]-2 for step in path if G.degree[step.frag_node]>2 ])

    def findBranchPath(self, subG:GraphFrags, source:FragNode) -> list[Pathstep]:
        if not self.canonize and not self.casual:
            return self.traverse_by_idx(list(subG._node), source)
        
        paths = self.find_paths(source, subG=subG)
        ## of course if it is just 1 ...
        if len(paths)==1:
            return paths[0]
        elif not self.canonize and self.casual:
            return secure_random.choice ( paths )
        
        paths=self.filter_paths(paths) # NOTE select only longest ones
        if len(paths)==1:
            return paths[0]
        
        paths = np.array(paths, dtype='O')
        ## Which path is the longest one ? Canonizing way ...
        numBranches=[self.getNumBranches(path, subG=subG) for path in paths]
        _maxNum=np.max(numBranches)
        if np.where( numBranches==_maxNum )[0].shape[0]==0 :
            return paths.tolist()[np.argmax(numBranches)]
        
        _minSum=self.getUniqueMinSum(paths)
        if not _minSum:## full symmetry !
            return paths.tolist()[0]
        else:
            return paths.tolist()[_minSum]
        
    def sortBranches(self, paths:list[list[Pathstep]]) -> list[int]:

        if len(paths) == 1:
            idxSorted=[0]
        elif not self.canonize and self.casual:
            idxSorted = list(range(len(paths)))
            secure_random.shuffle(idxSorted)
        elif not self.canonize and not self.casual:
            idxSorted = [i[0] for i in sorted(enumerate(paths), key=lambda x:x[1][0].frag_node.idx)]
        # else:
        # it's gonna retrieving for canonical sorting of branching paths

        sortedPaths=[]
        toSort=paths.copy()
        
        while toSort:
            numBranches=[self.getNumBranches(l) for l in toSort]
            _maxNum=np.max(numBranches)
            if np.where( numBranches==_maxNum )[0].shape[0]==1 :
                sortedPaths.append(toSort.pop(np.argmax(numBranches)))
                continue

            lenghts=np.array([len(l) for l in toSort])
            _maxNum=np.max(lenghts)
            if np.where( lenghts==_maxNum )[0].shape[0]==1 :
                sortedPaths.append(toSort.pop(np.argmax(lenghts) ))
                continue

            _minSum=self.getUniqueMinSum(toSort)
            if _minSum is not None:
                sortedPaths.append(toSort.pop(_minSum))
                continue

            sortedPaths+=toSort
            toSort=[]

        idxSorted = [paths.index(p) for p in sortedPaths]
        return idxSorted