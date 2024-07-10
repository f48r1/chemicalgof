import itertools
import numpy as np
import networkx as nx
from .tokens import TokenPath
import random
secure_random = random.SystemRandom()


def findLongestPath(G, source=None, lenMax=-1):
    terminalNodes = [node for node in G._node if G.degree[node]==1 and node!=source]
    # if the start of the path is not specified, all couples of terminal nodes are found
    if not source:
        terminalNodesCouple=[[(a,b),(b,a)] for a,b in itertools.combinations(terminalNodes,2)]
        terminalNodesCouple=itertools.chain(*terminalNodesCouple)
    # only couples of source and terminal nodes are found
    else:
        terminalNodesCouple=itertools.product([source],terminalNodes)

    ## path composed from 1 node ... :)
    if not terminalNodes:
        return np.array([[source]])

    paths=[]
    # to find the longest paths, still of unknown steps
    if lenMax<1:
        lenMax=0
        for a,b in terminalNodesCouple:
            path,*_=list(nx.all_simple_paths(G, a,b) )
            l=len(path)
            if l>lenMax:
                paths=[path]
                lenMax=l
            elif l==lenMax:
                paths.append(path)
    else: # to find only paths composed of n steps
        for a,b in terminalNodesCouple:
            path,*_=list(nx.all_simple_paths(G, a,b) )
            l=len(path)
            if l==lenMax:
                paths.append(path)

    return np.array(paths)

# path input is referred to the main path of reduced graph. So first and last element of list 'path' are terminal nodes
def findFirstBranchIdx(G, path):
    #we can ignore terminal nodes cause of they have just 1 deegree for sure.
    for i,node in enumerate(path[1:-1]):
        if G.degree[node]>2:
            return i+1
    return np.nan ## it is linear.

def getUniqueMinSum(paths):
    s=[0 for _ in paths ]
    for frags in zip(*paths):
        for idx, frag in enumerate(frags):
            s[idx] += frag.typeId
        minimum = min(s)
        _min = [pos for pos,element in enumerate(s) if element==minimum]
        if len(_min)==1 :
            return _min[0]
    return None ## No any path ranked as first

def getMultiMinSum(paths):
    s=np.zeros(len(paths), dtype=object )
    idxs=np.array(range(len(paths)))
    for frags in zip(*paths):
        IDs=np.array([f.typeId for f in frags], dtype=object)[idxs]
        s = [a+b  for a,b in zip (s,IDs)]
        minimum = min(s)
        _min = [pos for pos,element in enumerate(s) if element==minimum]
        if len(_min)==1 :
            return _min
        elif len(_min)<len(idxs):
            idxs = idxs[_min]
            s = [element for i,element in enumerate(s) if i in _min ]
    return idxs

### To delete maybe ###
def groupFragsByIdx(paths):
    sortedF=np.array([ sorted([frag.idx for frag in path]) for path in paths])
    groups=[ np.where( (sortedF==k).all(axis=1) )[0] 
            for k,_ in itertools.groupby(sortedF.tolist()) ]
    return groups

def groupFragsByType(paths):
    pathsFrags=[ tuple(str(f.typeId) for f in p) for p in paths ]
    groups=[ sorted(np.where( (np.array(pathsFrags)==k).all(axis=1))[0] ) 
            for k in set(pathsFrags) ]
    return groups


def buildBranches(mainDiG, mainG, mainPath, ascendent=None, canonizing=False):
    tokenP=TokenPath(mainPath)

    for pos,(n,t) in enumerate(zip(mainPath, tokenP)):
        if n.numPotAtomLinkers > 1:
            if pos==0 and ascendent!=None:
                dat=mainDiG.get_edge_data(n, ascendent)
                t.prec=str(dat["aB"]) if "stereo" not in dat else f"{dat['aB']}{dat['stereo']}"
            elif pos!=0 :
                dat=mainDiG.get_edge_data(n, mainPath[pos-1])
                t.prec=str(dat["aB"]) if "stereo" not in dat else f"{dat['aB']}{dat['stereo']}"

            if pos<len(mainPath)-1:
                dat=mainDiG.get_edge_data(n, mainPath[pos+1])
                t.succ=str(dat["aB"]) if "stereo" not in dat else f"{dat['aB']}{dat['stereo']}"

        # Secondary paths (branches)
        tmpNeighs=[neighbour for neighbour in mainG.neighbors(n) 
                   if neighbour not in mainPath and neighbour!=ascendent]
        if tmpNeighs:
            neighs=[]
            neighsT=[]
            for neigh in tmpNeighs:
                branch=getBranch(mainG, n, neigh)
                branchP=findBranchPath(mainG.subgraph(branch), neigh, canonizing)
                branchT=buildBranches(mainDiG, mainG, branchP, ascendent=n, canonizing=canonizing)
                neighs.append([x.frag for x in branchT])
                neighsT.append(branchT)
            
            if len(neighs)>1 and canonizing:
                idxSorted=sortBranches(mainG, neighs)
            elif len(neighs)>1 and not canonizing:
                idxSorted = list(range(len(neighs)))
                secure_random.shuffle(idxSorted)
            else:
                idxSorted=[0]
                
            for i in idxSorted:
                t.branches.append(neighsT[i])
                if n.numPotAtomLinkers > 1:
                    dat=mainDiG.get_edge_data(n, neighs[i][0])
                    t.branchesLinkers.append( str(dat["aB"]) if "stereo" not in dat else f"{dat['aB']}{dat['stereo']}" )

    return tokenP

    
def getBranch(G, core, source):
    subG=G.copy()
    subG.remove_edge(core, source)
    toDel=[x for x in nx.connected_components(subG) if source not in x]
    subG.remove_nodes_from(*toDel)
    return subG

def getNumBranches(G, path):
    return sum([G.degree[node]-2 for node in path if G.degree[node]>2 ])

def findBranchPath(subG, source, canonizing=False):
    longests=findLongestPath(subG, source)
    ## of course if it is just 1 ...
    if len(longests)==1:
        return np.array(*longests)
    elif not canonizing:
        return secure_random.choice ( longests )    
    
    ## Which path is the longest one ?
    numBranches=[getNumBranches(subG, path) for path in longests]
    _maxNum=np.max(numBranches)
    if np.where( numBranches==_maxNum )[0].shape[0]==0 :
        return longests[np.argmax(numBranches)]
    
    _minSum=getUniqueMinSum(longests)
    if not _minSum:## full symmetry !
        return longests[0]
    else:
        return longests[_minSum] 
    
def sortBranches(G, paths):
    sortedPaths=[]
    toSort=paths.copy()
    
    while toSort:
        numBranches=[getNumBranches(G, l) for l in toSort]
        _maxNum=np.max(numBranches)
        if np.where( numBranches==_maxNum )[0].shape[0]==1 :
            sortedPaths.append(toSort.pop(np.argmax(numBranches)))
        else:
            lenghts=np.array([len(l) for l in toSort])
            _maxNum=np.max(lenghts)
            if np.where( lenghts==_maxNum )[0].shape[0]==1 :
                sortedPaths.append(toSort.pop(np.argmax(lenghts) ))
            else:
                _minSum=getUniqueMinSum(toSort)
                if _minSum:
                    sortedPaths.append(toSort.pop(_minSum))
                else: ## indifferent
                    sortedPaths+=toSort
                    toSort=[]
                    
    return [paths.index(p) for p in sortedPaths]

    
# canonization of main path and secondary ones
def CanonicalGoF2Tokens(mainDiG):
    mainG = mainDiG.to_undirected()

    longests=findLongestPath(mainG)
        
    ## find index of first branching node for each path ##
    branchesIdxs=np.array([findFirstBranchIdx(mainG, p) for p in longests])


    ## we retain only path which node early branches
    if not np.isnan(branchesIdxs).all():
        _minBranch=np.nanmin(branchesIdxs)
        longests=longests[branchesIdxs==_minBranch]


    ## if only 1 path has early node branching
    if len(longests)==1:
        return buildBranches(mainDiG, mainG, longests[0], canonizing=True)


    branchesAmount = [ 
                        np.array([mainG.degree[branchingNode]-2
                        for branchingNode in path if mainG.degree[branchingNode]>2]) 
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

    grpByType=groupFragsByType(longests)

    if len(grpByType)>1:
        _minSum=getMultiMinSum(longests[ [i[0] for i in grpByType] ] )
        if _minSum is not None:
            _minSum=_minSum#.tolist()
            if np.iterable(_minSum):
                _inds=list ( itertools.chain( *[grpByType[ind] for ind in _minSum] ) )
            else:
                _inds=grpByType[_minSum]
            longests=longests[ _inds ]
            
            grpByType=groupFragsByType(longests)

    grpByType=groupFragsByType(longests)

    # indifferent
    if len(grpByType)==1 :
        return buildBranches(mainDiG, mainG, longests[0], canonizing=True)

    _minSum=getUniqueMinSum(longests[ [i[0] for i in grpByType] ] )

    if _minSum == None:
        ## should not happens ...
        print('Not predicted case !')
        return None

    return buildBranches(mainDiG, mainG, longests[_minSum], canonizing=True)

# return a random path. nNodesMainPath is setted deafult to the longest one
def GoF2Tokens(mainDiG, nNodesMainPath=-1) :
    mainG = mainDiG.to_undirected()

    if nNodesMainPath==-1:
        longests=findLongestPath(mainG)

    return buildBranches (mainDiG, mainG, secure_random.choice ( longests) )

# return more than 1 random path. it starts from the longests
def GoF2MoreTokens(mainDiG, nAugs=5):
    mainG = mainDiG.to_undirected()

    rets=[]
    strings=[]
    lenMax=-1
    while len(strings)<nAugs:
        longests=findLongestPath(mainG, lenMax=lenMax)
        if 1<=lenMax<3:
            break
        elif longests.size>0 and lenMax==-1:
            lenMax=len(longests[0])-1
        elif lenMax>1:
            lenMax-=1
        longests=longests.tolist()
        secure_random.shuffle ( longests )
        for path in longests:
            i=0
            canonizing = True
            while i<5:
                t=buildBranches(mainDiG, mainG, path, canonizing=canonizing)
                s=t.getString()
                if s not in strings:
                    strings.append(s)
                    rets.append(t)
                if len(strings)>=nAugs: break
                canonizing = False
                i+=1

    if len(rets)>nAugs:
        rets=rets[:nAugs]
    return rets
