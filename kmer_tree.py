from logging import addLevelName


class KMerTree ():
    def __init__(self, name):
        self.name = name
        self.children = {}
        self.reads = set()


    def add_kmer(self, kmer : str, read_name: str) -> None :
        if len(kmer) == 0 :
            self.reads.add(read_name)
            return None
        if kmer[0] not in self.children :
            self.children[kmer[0]] = KMerTree(kmer[0])
        self.children[kmer[0]].add_kmer(kmer[1:], read_name)

    def has_child(self, car : str) -> bool :
        return car in self.children

    def get_child(self, car : str) -> 'KMerTree' :
        if car not in self.children:
            raise KeyError("No child named ", car)
        return self.children[car]

    def get_last_node(self, kmer : str) -> 'KMerTree' or None :
        if len(kmer) == 0:
            return self
        if kmer[0] not in self.children :
            return None
        return self.children[kmer[0]].get_last_node(kmer[1:])

    def remove_kmer(self, kmer, read_name) :
        if len(kmer) == 0 :
            self.reads.remove(read_name)
            return None # on est dans le noeud qui porte le set
        target_child = self.children[kmer[0]]
        target_child.remove_kmer(kmer[1:], read_name)
        if not target_child.reads and not target_child.children :
            self.children.pop(kmer[0])
