import re
from .tokens import Token, TokenPath
from .utils import GetPotAtomLinkers
from rdkit import Chem
import numpy as np

# split to return frag followed by its binding atom idxs (if they are)
def splitByDots(s):
    s=list(s)
    i=0
    opened=0
    while i<len(s):
        if s[i]=="." and opened==0:
            s[i]=" "
        elif s[i]=="(" and (s[-1]=="." or s[-1]==">"):
            opened+=1
        elif s[i]==")" and s[-1]==".":
            opened-=1
        i+=1
    s="".join(s[:-1])
    return s.split(" ")

# split to return all branches of frag
def splitSecBranches(s):
    branches=[]
    while len(s)>1:
        i=0
        opened=0
        while opened==0:
            if s[i]=="(":
                opened+=1
            i+=1
        while True:
            if s[i]=="(" and (s[i-1]=="." or s[i-1]==">" or s[i-2:i]==".)"):
                opened+=1
            elif s[i]==")" and s[i-1]==".":
                opened-=1
                if opened==0:
                    branches.append(s[:i+1])
                    s=s[i+1:]
                    break
            i+=1
    return branches

# function employed by user to return tokenized reduced graph
def String2Tokens(s):
    mainP=splitByDots(s)
    frags=[]
    for f in mainP:
        if ".)" in f: ## branch
            branches=splitSecBranches(f)
            last=frags[-1]
            for b in branches:
                link=re.search(r"^<\d+[RS]?>",b)
                if link is not None:
                    b=b[link.end():]
                    last.branchesLinkers.append( link.group()[1:-1] )
                # else:
                #     last.branchesLinkers.append( None )

                b=b[1:-1]
                last.branches.append( String2Tokens(b) )
        else:
            pre=re.search(r"^<\d+[RS]?>",f)
            if pre:
                f=f[pre.end():]
                pre=pre.group()[1:-1]

            succ=re.search(r"<\d+[RS]?>$",f)
            if succ:
                f=f[:succ.start()]
                succ=succ.group()[1:-1]

            ## TODO
            ## canonicalize fragSMILES string to convert it to tokens and the directed graph    
            # if False:
            #     molDumm=prepareFragCanonization(f)
            #     smDumm=Chem.MolToSmiles(molDumm, canonize=False)
            #     canDumm=Chem.CanonizeSmiles(smDumm)
            #     if smDumm==canDumm:
            #         molCan, mapIdx = CanonizeFragWithDummies(molDumm)
            #         mapIdx[None]=None
            #         f=Chem.MolToSmiles(molCan)
            #         pre=mapIdx[pre]
            #         succ=mapIdx[succ]

            #frag=FragNode.fromSmiles( f )
            fragT=Token(f)
            fragT.prec=pre
            fragT.succ=succ
            frags.append(fragT)

    return TokenPath(frags)

def Sequence2String(arr):
    ids=set()
    i=0
    while i<len(arr):
        if not any (arr[i].startswith(_) for _ in ["<","(",")"]): # frammento
            if len(GetPotAtomLinkers(arr[i]))>1:
                if i>0:
                    if arr[i-2]!="(":
                        ids.add(i-1)
                if i<len(arr)-1:
                    ids.add(i+2)
            else:
                ids.add(i+1)
        elif arr[i]==")":
            ids.add(i)
            if i<len(arr)-1:
                if arr[i+1]!="(" and ( not arr[i+1].startswith("<") and arr[i+1]!="("):
                    ids.add(i+1)
        i+=1

    ids.add(len(arr))
    ret= np.insert(arr, list(ids), ".")
    return "".join(ret)

def String2Sequence(s):
    return String2Tokens(s).getSequence()

def Sequence2Smiles(x):
    try:
        x=Sequence2String(x)
        x=String2Tokens(x)
        x=x.getGraph().getMol()
        return Chem.MolToSmiles(x)
    # except Exception as e:
    #     # print(e)
    except:
        return None