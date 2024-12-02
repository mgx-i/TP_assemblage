from kmer_tree import *
import pytest

def test_kmer_initialisation () :
    tree = KMerTree("root")
    tree.add_kmer("ACGT")
    assert tree.has_child("A") is True
    assert tree.has_child("G") is False
    anode1 = tree.get_child("A")
    cnode2 = anode1.get_child("C")
    gnode3 = cnode2.get_child("G")
    tnode4 = gnode3.get_child("T")