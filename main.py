from kmer_tree import *

def reverse_sequence (s : str) -> str :
    reversed_nucleotide = {"A": "T",
                     "T": "A",
                     "G": "C",
                     "C": "G"}
    return "".join(reversed([reversed_nucleotide[nucleotide] for nucleotide in s]))


def read_fastq():
    #version test
    reversed_nucleotide = {"A": "T",
                     "T": "A",
                     "G": "C",
                     "C": "G"}
    l_ = [
        "AAATGCGATCCGATAGGCA",
        "GATAGGCAAGGTGAGCTAG",
        "GCTAGGGGATTCACTGATT",
        "ACTGATTAAGGGTCGACGA",
        "AGGGTCACCCGATCGATCA",
        "CGATCGATGCGGGAGTTTC"
    ]

    for i, reads in enumerate(l_):
        yield reads, f"{i}f"
        yield reverse_sequence(reads), f"{i}r"


def fill_kmer_tree (tree : KMerTree, read : str, tag : str, k : int) :
    for i in range(len(read) - k - 1) :
        tree.add_kmer(read[i:i+k], tag)


def pick_kmer (tree : KMerTree) -> str :
    #choisir le premier kmer venu
    kmer = ''
    node = tree
    while node.children :
        nucleotide = list(node.children.keys())[0]
        kmer += nucleotide
        node = node.children[nucleotide] # va pas du tout
    return kmer


def extend_kmer (contig : str, k : int, tree : KMerTree) -> str :
    end_kmer = contig[len(contig) - k + 1:]
    last_node = tree.get_last_node(end_kmer)

    while last_node :
        if len(last_node.children) == 1:
            new_nucleotide = list(last_node.children.keys())[0]
            contig += new_nucleotide

            forward_kmer = end_kmer + new_nucleotide
            reverse_kmer = reverse_sequence(forward_kmer)
            tree.remove_kmer(forward_kmer)
            tree.remove_kmer(reverse_kmer)
        else :
            return contig

        end_kmer = contig[len(contig) - k + 1:]
        last_node = tree.get_last_node(end_kmer)

    return contig


def assemble (k : int) -> None :
    kmer_tree = KMerTree(name="root")
    for read, tag in read_fastq():
        fill_kmer_tree(kmer_tree, read, tag, k)

    f = open("contigs.txt", "w")
    i = 1
    while kmer_tree.children : # tant qu'il y a des kmers dans le tree
        first_kmer = pick_kmer(kmer_tree)

        # étend le kmer vers la droite
        intermediate_contig = extend_kmer(first_kmer, k, kmer_tree)
        # on retourne le contig pour l'étendre vers la gauche :
        reverse_intermediate_contig = reverse_sequence(intermediate_contig)
        contig = extend_kmer(reverse_intermediate_contig, k, kmer_tree)

        f.write(f">Contig {i}\n")
        f.write(contig)
        f.write('\n')
        i += 1

    f.close()

assemble(5)
