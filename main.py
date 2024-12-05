from kmer_tree import *
import timeit

def reverse_sequence (s : str) -> str :
    """Give the reverse complement of DNA sequence.
     :param str s: A nucleotide sequence
     :return str: The reverse complement of s (5'-3').
    """
    reversed_nucleotide = {"A": "T",
                     "T": "A",
                     "G": "C",
                     "C": "G"}
    return "".join(reversed([reversed_nucleotide[nucleotide] for nucleotide in s]))


def read_fastq(file_name: str) :
    """Retrieve reads from a fastQ.
    :param str file_name: A path that lead to a fastQ file.
    :return: Python generator which contain all fastQ's sequences and a identifier for the read:
        -> (Read1, 1)
    """
    f = open(file_name, "r")
    next_is_fastq = True
    i=0
    line = f.readline()

    while line:
        line = f.readline()
        if next_is_fastq:
            yield line, str(i)
            i+=1
        if line and line[0] == '@':
            next_is_fastq = True
        else:
            next_is_fastq = False
    f.close()

def fill_kmer_tree (tree : KMerTree, read : str, tag : str, k : int) :
    """Add all kmer contained in a read inside a tree.
    :param KMerTree tree: A tree of kmer.
    :param str read: A Read from a fastQ file. It's length should be greater or equal to k.
    :param str tag: A tag for those kmer. Usually the origin of this read.
    :param int k : An odd number greater than 0. Kmer size.
    """
    for i in range(len(read) - k - 1) :
        tree.add_kmer(read[i:i + k], tag + "s")
        tree.add_kmer(reverse_sequence(read[i:i + k]), tag + "r")


def pick_kmer (tree : KMerTree) -> str :
    """Return the first kmer added to a KMerTree
    :param KMerTree tree: A tree of kmer.
    :return str: A kmer from KMerTree.
    """
    # Select the first kmer added to the tree.
    kmer = ''
    node = tree

    while node.children : # If the node have no children, the node is a leaf.
        # Reconstruct kmer's sequence.
        nucleotide = list(node.children.keys())[0]
        kmer += nucleotide
        node = node.children[nucleotide] # va pas du tout
    return kmer


def extend_kmer (contig : str, k : int, tree : KMerTree) -> str :
    """Extend a contig using a KMerTree.
    :param str contig: A contig that should be extended.
    :param int k : An odd number greater than 0. Kmer size.
    :param KMerTree tree: A tree of kmer
    :return str: A contig.
    """

    # Set up the loop
    end_kmer = contig[len(contig) - k + 1:] # Suffix of "contig" sized k-1
    last_node = tree.get_last_node(end_kmer) # Node that represent the suffix of "contig" sized k-1

    while last_node :
        if len(last_node.children) == 1:
            # Elongate the contig
            new_nucleotide = list(last_node.children.keys())[0]
            contig += new_nucleotide


            # Delete used k-mer to avoid infinite while loop
            forward_kmer = end_kmer + new_nucleotide
            reverse_kmer = reverse_sequence(forward_kmer)
            tree.remove_kmer(forward_kmer)
            # Since we add the reverse complete along with the standard string, we have
            # to remove it.
            tree.remove_kmer(reverse_kmer)

        else :
            # Either no or multiple children.
            return contig

        # Contunue the loop
        end_kmer = contig[len(contig) - k + 1:]
        last_node = tree.get_last_node(end_kmer)

    return contig


def assemble (path: str, k : int) -> int :
    """Realise an assembly using reads inside a uniq fastQ file.
    :param str path: A path that lead to a fastQ file.
    :param int k : An odd number greater than 0. Kmer size.
    All contig are wrote inside contigs_{k}.fa where {k} is str(k)
    """
    # Fill the kmer tree
    kmer_tree = KMerTree(name="root")
    for read, tag in read_fastq(path):
        fill_kmer_tree(kmer_tree, read, tag, k)

    # Prepare the loop
    f = open(f"contigs_{k}.fa", "w")
    i = 1

    while kmer_tree.children : # While there is at least one kmer in the tree.
        # Select a kmer in order to start the process.
        first_kmer = pick_kmer(kmer_tree)
        kmer_tree.remove_kmer(first_kmer)
        kmer_tree.remove_kmer(reverse_sequence(first_kmer))

        # contig elongation in the 5'-3' direction
        intermediate_contig = extend_kmer(first_kmer, k, kmer_tree)
        # contig elongation in the 3'-5' direction
        reverse_intermediate_contig = reverse_sequence(intermediate_contig)
        contig = extend_kmer(reverse_intermediate_contig, k, kmer_tree)

        # Output the contig.
        f.write(f">{i}__len__{len(contig)}\n")
        f.write(contig)
        f.write('\n')
        i += 1

    f.close()
    return i


def contigs_for_k(path, kn) :
    """Internal test function. Test all possible odd k-mer size from 3 to kn using the reads in a file "path"."""
    f = open("contigs_stats.txt", "w")
    f.write("k \t nb contigs\n")
    for k in range(3,kn+1,2):
        n_contigs = assemble(path, k)
        f.write(f"{k}\t{n_contigs}\n")
    f.close()


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Two arguments needed : \n\t- fastq file path \n\t- kmer size")
        exit(1)

    try:
        print(f"Assembly from file {sys.argv[1]} with k={sys.argv[2]}.")
        print(f"Executed in {timeit.timeit(lambda : assemble(path=sys.argv[1], k=int(sys.argv[2])), number=1)}s.")

    except Exception as E:
        print(f"Unexpected exception : {E}")
        exit(2)
