from .classes import FragNode, GraphFrags, DiGraphFrags
import networkx as nx, itertools

# classes for tokenized graph handling

class Token(str):       
    def __new__(cls, frag):
        
        if type(frag)==str:
            frag=FragNode.fromSmiles(frag)
            
        if frag.suffStereo is not None:
            instance = super().__new__(cls, f"{frag.smiles}|{frag.suffStereo}")
        else:
            instance = super().__new__(cls, f"{frag.smiles}")
        instance.frag=frag
            
        instance.branches=[]
        instance.prec=None
        instance.succ=None
        instance.branchesLinkers=[]
        return instance
    
    def __repr__(self):
        b=f"({len(self.branches)})" if self.branches else ""
        return self+b
    
# tokenPath can be just one token node, but it is still a path !
class TokenPath(list):
   
    def __init__(self, iterable):
        super().__init__(Token(item) if type(item)!=Token else item for item in iterable)

    def getSequence(self):
        ret=[]
        for t in self:
            if t.prec!=None:
                ret.append(f"<{t.prec}>")
            ret.append(str(t))
            if t.succ!=None:
                ret.append(f"<{t.succ}>")
            if t.branches:
                if t.branchesLinkers:
                    for p,l in zip(t.branches, t.branchesLinkers):
                        ret.append(f"<{l}>")
                        ret+=["(",*p.getSequence(),")"]
                else:
                    for p in t.branches:
                        ret+=["(",*p.getSequence(),")"]
        return ret

    def getString(self):
        ret=""
        for t in self:
            if t.prec is not None:
                ret+=f"<{t.prec}>"
            ret+=str(t)
            if t.succ is not None:
                ret+=f"<{t.succ}>"
            ret+="."
            if t.branches:
                if t.branchesLinkers:
                    for p,l in zip(t.branches, t.branchesLinkers):
                        ret+=f"<{l}>"+"("+p.getString()+")"
                else:
                    for p in t.branches:
                        ret+="("+p.getString()+")"
                ret+="."
        return ret
    
    def getGraph(self):
        diG=DiGraphFrags()
        nodes=[n.frag for n in self]
        diG.add_nodes_from(nodes)

        for (nPre,tPre),(nSucc, tSucc) in itertools.pairwise(zip(nodes,self)):
            kwargs={}
            if tPre.succ!=None:
                if tPre.succ.endswith("S") or tPre.succ.endswith("R"):
                    kwargs["aB"]=int(tPre.succ[:-1])
                    kwargs["stereo"]=tPre.succ[-1]
                else:
                    kwargs["aB"]=int(tPre.succ)
            else:
                kwargs["aB"]=nPre.PotAtomLinkers[0]
                if nPre.suffStereo!=None:
                    kwargs["stereo"]=nPre.suffStereo
                
            diG.add_edge(nPre, nSucc, **kwargs)
            
            kwargs={}
            if tSucc.prec!=None:
                if tSucc.prec.endswith("S") or tSucc.prec.endswith("R"):
                    kwargs["aB"]=int(tSucc.prec[:-1])
                    kwargs["stereo"]=tSucc.prec[-1]
                else:
                    kwargs["aB"]=int(tSucc.prec)
            else:
                kwargs["aB"]=nSucc.PotAtomLinkers[0]
                if nSucc.suffStereo!=None:
                    kwargs["stereo"]=nSucc.suffStereo
                
            diG.add_edge(nSucc, nPre, **kwargs)
        
        ## Branches of each node
        for node,tok in zip(nodes,self):
            
            for i,p in enumerate(tok.branches):
                H=p.getGraph()
                diG=nx.compose( diG, H )
                other=p[0]
                
                kwargs={}
                if tok.branchesLinkers:
                    if tok.branchesLinkers[i].endswith("S") or tok.branchesLinkers[i].endswith("R"):
                        kwargs["aB"]=int(tok.branchesLinkers[i][:-1])
                        kwargs["stereo"]=tok.branchesLinkers[i][-1]
                    else:
                         kwargs["aB"]=int(tok.branchesLinkers[i])
                else:
                    if node.suffStereo!=None:
                        kwargs["stereo"]=node.suffStereo
                    kwargs["aB"]=node.PotAtomLinkers[0]
                    
                diG.add_edge(node, other.frag, **kwargs)
                
                
                kwargs={}
                if other.prec:
                    if other.prec.endswith("S") or other.prec.endswith("R"):
                        kwargs["aB"]=int(other.prec[:-1])
                        kwargs["stereo"]=other.prec[-1]
                    else:
                        kwargs["aB"]=int(other.prec)
                else:
                    if other.frag.suffStereo!=None:
                        kwargs["stereo"]=other.frag.suffStereo
                    kwargs["aB"]=other.frag.PotAtomLinkers[0]
                diG.add_edge(other.frag, node, **kwargs)
                
                ### update parent graph for each node
                for _ in diG._node: _.setParent(diG)
                
        return diG