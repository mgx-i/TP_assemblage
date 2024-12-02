from kmer_tree import *
import pytest

def create_tree () :
    tree = KMerTree("root")
    tree.add_kmer("ACGT", "1s")
    tree.add_kmer("GCTT", "1s")
    tree.add_kmer("ACTA", "2s")
    tree.add_kmer("ACGT", "2a")
    tree.add_kmer("TTTT", "2a")
    return tree


def test_kmer_initialisation () :
    tree = KMerTree("root")

    tree.add_kmer("ACGT", "1a")
    tree.add_kmer("ACGT", "2s")
    assert tree.has_child("A") is True
    assert tree.has_child("G") is False
    anode1 = tree.get_child("A")
    cnode2 = anode1.get_child("C")
    gnode3 = cnode2.get_child("G")
    tnode4 = gnode3.get_child("T")
    assert tnode4.reads == {"1a", "2s"}

    tree.add_kmer("ACTA", "1a")
    tnode3 = cnode2.get_child("T")
    anode4 = tnode3.get_child("A")
    assert anode4.reads == {"1a"}


def test_get_last_node () :
    tree = create_tree()
    anode1 = tree.get_child("A")
    cnode2 = anode1.get_child("C")
    gnode3 = cnode2.get_child("G")
    tnode4 = gnode3.get_child("T")
    assert tree.get_last_node("ACGT") == tnode4
    assert tree.get_last_node("ACG") == gnode3


