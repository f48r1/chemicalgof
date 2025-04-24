from .traverse import Pathstep
from .gof import FragNode, DiGraphFrags
from .traverse import Traverser, secure_random # NOTE this is needed for the random choise of the main path.

import numpy as np
import itertools

class Writer:
    def __init__(self, DiG:DiGraphFrags):
        self.DiG = DiG

    def write_bond(self, source: FragNode, target: FragNode):
        if source.numPotAtomLinkers <=1:
            return ''
        
        data_bond = self.DiG.get_edge_data(source, target)

        # if data_bond['aB'] is None: # [ ] set None for the 'aB' node with 1 only linker ? Next versions maybe
            # return ''

        bond_str = str(data_bond['aB'])
        if data_bond['stereo'] is not None : bond_str+=str(data_bond['stereo'])

        return f'<{bond_str}>'

    def write_fragsmiles(self, path:list[Pathstep], ascendent_node: FragNode | None = None, sep='.') -> str:

        string = ''
        if ascendent_node is not None:
            string += self.write_bond(path[0].frag_node, ascendent_node)

        # NOTE since first node cant branching, its fragSMILES is written.
        # NOTE Same result if last node was written out of the cycle for
        string += str(path[0].frag_node)

        for prec,succ in itertools.pairwise(path):
            prec_node = prec.frag_node; succ_node = succ.frag_node
            string+=self.write_bond(prec_node, succ_node) + sep

            if prec.branches:
                for branching in prec.branches:
                    string+=self.write_bond(prec_node, branching[0].frag_node) + '('
                    # [ ] do we need seriusly to add sep before closing bracket ?? It is just because of re.findall (tokenization rule)
                    string+=self.write_fragsmiles(branching, ascendent_node=prec_node) + sep + ')' + sep

            string+=self.write_bond(succ_node, prec_node)
            string += str(succ_node)

        return string
    
# canonization of main path and secondary ones
def CanonicalGoF2fragSMILES(DiG:DiGraphFrags) -> str | None: # NOTE None is for exceptional cases of bugs !
    traverser = Traverser(DiG, canonize=True, casual=False)
    writer = Writer(DiG)

    if DiG.number_of_nodes() == 1:
        pathway=[Pathstep(node) for node in DiG._node]
        return writer.write_fragsmiles(pathway)

    longests=traverser.find_paths()
    longests=traverser.filter_paths(longests)

    if longests is None:
        print('Latest BUG ! SMARTS pattern doesn\'t recognize ciclyc structure and then graph is cyclic !' )
        return None
    
    longests = np.array(longests, dtype='O') # Convert into array to more easy handling
        
    ## find index of first branching node for each path ##
    branchesIdxs=np.array([traverser.findFirstBranchIdx(path) for path in longests])

    ## we retain only path which node early branches
    if not np.isnan(branchesIdxs).all():
        _minBranch=np.nanmin(branchesIdxs)
        longests=longests[branchesIdxs==_minBranch]

    ## if only 1 path has early node branching
    if len(longests)==1:
        canonical_path = traverser.buildBranches(longests[0])

        return writer.write_fragsmiles(canonical_path)


    branchesAmount = [ 
                        np.array([traverser.G.degree[step.frag_node]-2
                        for step in path if traverser.G.degree[step.frag_node]>2]) 
                        for path in longests
                    ]

    _lens=[len(x) for x in branchesAmount]
    if len(set(_lens))>1:
        iMax=np.max(_lens)
        for pos,l in enumerate(_lens):
            if l<iMax:
                branchesAmount[pos]=np.pad(branchesAmount[pos], (0, iMax-l) )
        
    branchesAmount=np.array( branchesAmount , dtype=np.int64)
    
    for col in branchesAmount.T:
        _maxBr=np.max(col)
        if np.argwhere( col==_maxBr ).size != len(longests):
            longests=longests[col==_maxBr]
            break

    grpByType=traverser.groupFragsByType(longests)

    if len(grpByType)>1:
        _minSum=traverser.getMultiMinSum(longests[ [i[0] for i in grpByType] ] )
        if _minSum is not None:
            _minSum=_minSum#.tolist()
            if np.iterable(_minSum):
                _inds=list ( itertools.chain( *[grpByType[ind] for ind in _minSum] ) )
            else:
                _inds=grpByType[_minSum]
            longests=longests[ _inds ]
            
            grpByType=traverser.groupFragsByType(longests)

    grpByType=traverser.groupFragsByType(longests)

    # indifferent
    if len(grpByType)==1 :
        canonical_path = traverser.buildBranches(longests.tolist()[0])
        return writer.write_fragsmiles(canonical_path)

    _minSum=traverser.getUniqueMinSum(longests[ [i[0] for i in grpByType] ] )

    if _minSum == None:
        ## should not happens ...
        print('Not predicted case !')
        return None

    canonical_path = traverser.buildBranches(longests[_minSum])

    return writer.write_fragsmiles(canonical_path)

def RandomGoF2FragSMILES(DiG:DiGraphFrags):
    traverser = Traverser(DiG, canonize=False, casual=True)
    writer = Writer(DiG)

    paths = traverser.find_paths()
    main_path = secure_random.choice ( paths )
    main_path = traverser.buildBranches(main_path)

    return writer.write_fragsmiles(main_path)


def OrderedGoF2fragSMILES(DiG:DiGraphFrags):
    traverser = Traverser(DiG, canonize=False, casual=False)
    writer = Writer(DiG)

    main_path = traverser.traverse_by_idx()
    main_path = traverser.buildBranches(main_path)

    return writer.write_fragsmiles(main_path)

def GoF2fragSMILES(DiG:DiGraphFrags, canonize=False, random=False) -> str:
    """Convert reduced graph into fragSMILES representation. Graph is traversed in a canonical way, randomly or ordered (if canonize and random are both False).

    Args:
        DiG (DiGraphFrags): _description_
        canonize (bool, optional): Traversing reduced graph in a canonical way. Defaults to False.
        random (bool, optional): Traversing reduced graph randomly. Defaults to False.

    Returns:
        str: fragSMILES representation.
    """
    if canonize:
        return CanonicalGoF2fragSMILES(DiG)
    elif random:
        return RandomGoF2FragSMILES(DiG)
    else:
        return OrderedGoF2fragSMILES(DiG)