# [ ] A few exceptions could be developed to replace ValueError in other files.

class ConnectorIndex(Exception):
     
    def __init__(self, fragsmiles:str, index:int):
        self.fragsmiles = fragsmiles
        self.index = index
        self.message = f'Bond Error : Connecting atom index not specified for fragment: {self.fragsmiles} (index element #{self.index}).'
        super().__init__(self.message)

class InvalidChirality(Exception):
     
    def __init__(self, atom_symbol:str | None = None, index:int | None = None, expected:str | None = None):
        self.atom_symbol = atom_symbol
        self.index = index
        self.expected = expected
        if atom_symbol is not None  and index is not None and expected is not None:
            self.message = f'Chirality Error : Expected "{self.expected}" chirality label for atom {self.atom_symbol} #{self.index} is not allowed.'
        else:
            self.message = f'Chirality Error : Expected chirality label is not allowed.'
        super().__init__(self.message)