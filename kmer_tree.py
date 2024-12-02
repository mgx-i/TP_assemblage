from logging import addLevelName


class KMerTree ():
    def __init__(self, name):
        self.name = name
        self.children = {}


    def add_kmer(self, kmer : str) -> None :
        if len(kmer) == 0 :
            return None
        if kmer[0] not in self.children :
            self.children[kmer[0]] = KMerTree(kmer[0])
        self.children[kmer[0]].add_kmer(kmer[1:])

    def has_child(self, car : str) -> bool :
        return car in self.children

    def get_child(self, car : str) -> 'KMerTree' :
        if car not in self.children:
            raise KeyError("No child named ", car)
        return self.children[car]