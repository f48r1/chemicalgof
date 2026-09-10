from .gof import DiGraphFrags, FragNode, Fragment
from .exceptions import ConnectorIndex
from .utils import GetPotAtomLinkers
import re

_compiler_split = re.compile(r'(?<=>|.)\(|\)(?=.)|<[0-9]+[RSrs]?>|(?<=>|.)[^\.<>]+(?=<|.|\).)')

_compiler_split_atom_level = re.compile(
    r"\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|\*|\$|\%[0-9]{2}|(?<!<)[0-9](?!>)|<[0-9]+[A-Z]*>|[0-9]+[RS]|\|"
)
def split_fragsmiles(fragsmiles:str, atom_level=False) -> list[str]:
    """split fragSMILES string represention into list of tokens.

    Args:
        fragsmiles (str): fragSMILES representation

    Returns:
        list[str]: fragSMILES tokenized into list of strings.
    """
    if not atom_level:
        string = '<b>'+fragsmiles+'<e>'
        splitted = _compiler_split.findall(string)[1:-1]
    else:
        splitted = _compiler_split_atom_level.findall(fragsmiles)

    return splitted


class Parser:

    def __init__(self):
        self.DiG = DiGraphFrags()

    @staticmethod
    def setup_fragment(fragment_element:str) -> Fragment:
        # TODO
        if '|' in fragment_element:
            smiles,suffix = fragment_element.split('|')
            PotAtomLinkers = GetPotAtomLinkers(smiles)
            if len(PotAtomLinkers) == 1:
                chirality={PotAtomLinkers[0]:suffix}
            else:
                matches = re.findall(r'([0-9]+)([RSrs]{1})',suffix)
                chirality = {int(k):v for k,v in matches}
        else:
            smiles = fragment_element
            chirality = {}
        
        return Fragment(smiles=smiles, chirality=chirality)
    
    @staticmethod
    def setup_edge(edge_element:str) -> dict:
        matches = re.match(r'<(?P<aB>[0-9]+)(?P<stereo>[RSrs]?)',edge_element)
        attributes = matches.groupdict()
        if not attributes['stereo']:
            attributes['stereo'] = None
        attributes['aB'] = int(attributes['aB'])
        return attributes
    
    @staticmethod
    def extract_branching(sequence:list[str], starting_idx:int) -> tuple[list[str], list[str]]:
        main_sequence = sequence[:starting_idx]
        branching_sequence = sequence[starting_idx+1:]

        opened = 1
        closed = 0
        i=0
        while opened != closed and i<len(branching_sequence):
            element = branching_sequence[i]
            if element == '(':
                opened+=1
            elif element == ')':
                closed+=1
            i+=1
        
        if opened != closed:
            raise ValueError('Branching Error : Some branching path has not been closed or opened')

        main_sequence = main_sequence + branching_sequence[i:]
        branching_sequence = branching_sequence[:i-1]

        return main_sequence, branching_sequence

    def parse(self, sequence:list[str], ascendent_node:FragNode | None = None, ascendent_edge:dict | None = None) -> DiGraphFrags:

        if ascendent_edge is None:
            ascendent_edge = {}
        discendent_edge = {}

        idx_token_inline=0
        while idx_token_inline<len(sequence):
            element = sequence[idx_token_inline]

            if element == '(' and (ascendent_node is None or ascendent_edge is None):
                raise ConnectorIndex(fragsmiles='Undefined', index=idx_token_inline-1)
            
            elif element == '(':
                sequence, branching = self.extract_branching(sequence, idx_token_inline)
                if not discendent_edge and not ascendent_edge:
                    raise ConnectorIndex(fragsmiles=ascendent_node.smiles, index=idx_token_inline-1)

                # NOTE discendent_edge empty means hat branching node has not specified connecting atom
                branching_edge = {**discendent_edge} if discendent_edge else {**ascendent_edge}
                _ = self.parse(branching, ascendent_node=ascendent_node, ascendent_edge=branching_edge)

                discendent_edge.clear()
                idx_token_inline-=1 # NOTE branching is replace by next element

            elif not any (element.startswith(_) for _ in ["<","(",")"]): # fragment
                fragment = self.setup_fragment(element)
                node = FragNode(fragment)
                self.DiG.add_node(node)

                if ascendent_node is not None:
                    if ascendent_edge and discendent_edge and fragment.num_connector == 1: # multi atom bonding
                        raise ValueError('Bond Error : single atom bonding has also specificied connector index')
                    
                    elif ascendent_edge and not discendent_edge and fragment.num_connector > 1: # single atom bonding
                        raise ConnectorIndex(fragsmiles=fragment.smiles, index=idx_token_inline)
                    
                    elif not ascendent_edge:
                        raise ConnectorIndex(fragsmiles=fragment.smiles, index=idx_token_inline)
                    
                    elif not discendent_edge and fragment.num_connector == 1:
                        discendent_edge['aB'] = list(fragment.connector_idxs.keys())[0]
                        discendent_edge['stereo'] = None
                    
                    self.DiG.add_edge(node, ascendent_node, **discendent_edge)
                    self.DiG.add_edge(ascendent_node, node, **ascendent_edge)

                    discendent_edge.clear()
                    ascendent_edge.clear()

                if fragment.num_connector == 1:
                    ascendent_edge['aB'] = list(fragment.connector_idxs.keys())[0]
                    ascendent_edge['stereo'] = None

                ascendent_node = node

            elif element.startswith('<') and element.endswith('>'): # bond
                edge_attr = self.setup_edge(element)

                if ascendent_edge and discendent_edge:
                    raise ValueError('Bond Error: many consecutive connecting atoms')
                
                elif ascendent_edge and not discendent_edge:
                    discendent_edge.update(edge_attr)
                elif not ascendent_edge:
                    ascendent_edge.update(edge_attr)

            elif element == ')':
                raise ValueError('Branching Error : Closing bracket for no any branching path')

            idx_token_inline+=1
        return self.DiG
    
def Sequence2GoF(sequence:list[str]) -> DiGraphFrags:
    obj = Parser()
    DiG = obj.parse(sequence)

    return DiG
    
def fragSMILES2GoF(fragsmiles:str) -> DiGraphFrags:

    sequence = split_fragsmiles(fragsmiles)

    return Sequence2GoF(sequence)