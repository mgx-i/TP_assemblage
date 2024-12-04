from sys import displayhook

from kmer_tree import *

def reverse_sequence (s : str) -> str :
    reversed_nucleotide = {"A": "T",
                     "T": "A",
                     "G": "C",
                     "C": "G"}
    return "".join(reversed([reversed_nucleotide[nucleotide] for nucleotide in s]))


def read_fastq(file_name) :
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
    for i in range(len(read) - k - 1) :
        tree.add_kmer(read[i:i + k], tag + "s")
        tree.add_kmer(reverse_sequence(read[i:i + k]), tag + "r")


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
            return contig, "Multiple_possibilities"

        end_kmer = contig[len(contig) - k + 1:]
        last_node = tree.get_last_node(end_kmer)

    return contig, "No_matching_kmer"


def assemble (path: str, k : int) -> int :
    kmer_tree = KMerTree(name="root")
    for read, tag in read_fastq(path):
        fill_kmer_tree(kmer_tree, read, tag, k)

    f = open(f"contigs_{k}.txt", "w")
    i = 1
    while kmer_tree.children : # tant qu'il y a des kmers dans le tree
        first_kmer = pick_kmer(kmer_tree)
        kmer_tree.remove_kmer(first_kmer)
        kmer_tree.remove_kmer(reverse_sequence(first_kmer))

        # étend le kmer vers la droite
        intermediate_contig, right_stop = extend_kmer(first_kmer, k, kmer_tree)
        # on retourne le contig pour l'étendre vers la gauche :
        reverse_intermediate_contig = reverse_sequence(intermediate_contig)
        contig, left_stop = extend_kmer(reverse_intermediate_contig, k, kmer_tree)

        f.write(f">Contig {i} {right_stop}-{left_stop}\n")
        f.write(contig)
        f.write('\n\n')
        i += 1

    f.close()
    return i


def contigs_for_k(path, kn) :
    f = open("contigs_stats.txt", "w")
    f.write("k \t nb contigs\n")
    for k in range(3,kn+1,2):
        n_contigs = assemble(path, k)
        f.write(f"{k}\t{n_contigs}\n")
    f.close()


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Two arguments needed : \n\t- fastq file path \n\t- kmer size")
        exit(1)

    try:
        assemble(path=sys.argv[1], k=int(sys.argv[2]))
    except Exception as E:
        print(f"Unexpected exception : {E}")
        exit(2)
