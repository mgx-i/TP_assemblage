from logging import addLevelName


class KMerTree ():
    """
    Represent a number of kmer using a tree structure.
    Same kmers from different reads are merged together but their origin are kept in  self.reads. When a kmer is
    present multiple time inside an uniq read, this kmer isn't stored.

    self.name: str -> The nucleotide represented by this node.
    self.children: dict[str: KMerTree] -> Contain node's child. Child's names usually are a nucleotide.
    self.reads: set[str] -> A set that remember in which reads a kmer is from.
    """

    def __init__(self, name):
        """
        Initialise KMerTree.
        :param str name: A name for this KMerTree
        """
        self.name: str = name
        self.children: dict[str:KMerTree] = {}
        self.reads: set[str] = set()


    def add_kmer(self, kmer : str, read_name: str) -> None :
        """
        Add a kmer inside a KMerTree. Please only use this command on the KMerTree's root.
        
        Be aware that all kmer in a KMerTree should have the same length.

        :param str kmer: A kmer
        :param str read_name: A string that identify the origin of this kmer
        """
        if len(kmer) == 0 :
            # We reached the end of this kmer
            self.reads.add(read_name)   # Save kmer identifier in the leaf.
            return None
        if kmer[0] not in self.children :
            # Current node do not have this nucleotide as child. We have to generate a nex node
            # => We should add a new branch to the tree.
            self.children[kmer[0]] = KMerTree(kmer[0])
        self.children[kmer[0]].add_kmer(kmer[1:], read_name)   # Continue the insertion of this kmer

    def has_child(self, car : str) -> bool :
        """
        Does this node have a specific child.
        :param str car: the name of the child that should be child of this node.
        :return bool: Says whether a node is a child of the current node.
        """
        return car in self.children

    def get_child(self, car : str) -> 'KMerTree' :
        """
        Return a child belonging to the current node using its name.
        :param str car: Name of a child belonging to the current node.
        :return KMerTree: The targeted child.
        :raise KeyError: When the 'car' does not match any child name.
        """
        if car not in self.children:
            raise KeyError("No child named ", car)
        return self.children[car]

    def get_last_node(self, kmer : str) -> 'KMerTree' or None :
        """
        Return the leaf that correspond to a kmer.

        :param str kmer: A kmer.
        :return KMerTree: None when kmer is not in the KMerTree. KMerTree's node otherwise.
        """
        if len(kmer) == 0:
            return self
        if kmer[0] not in self.children :
            return None
        return self.children[kmer[0]].get_last_node(kmer[1:])

    def remove_kmer(self, kmer, read_name) :
        """
        Remove a kmer from the tree.
        A kmer is completely removed from the tree when all kmer's origin are deleted. (If a kmer comes from
        'read1' and 'read2', you should use remove_kmer(kmer, 'read1') and remove_kmer(kmer, 'read2').
        :param str kmer: A kmer in the tree
        :param read_name: A tag related to kmer's origin (read tag)
        """
        if len(kmer) == 0 :
            # We are in the node that contain the origin of this kmer
            self.reads.remove(read_name)
            return None
        target_child = self.children[kmer[0]]
        target_child.remove_kmer(kmer[1:], read_name)
        if not target_child.reads and not target_child.children :
            # This child have does not contain any valuable data. It must be exterminated.
            self.children.pop(kmer[0])
