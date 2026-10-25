import pytest
from genome import Genome

@pytest.fixture
def def_genome():
    """default genome for basic case"""
    return Genome("genome1", "ATCGATCGA")


def test_add_kmer(def_genome):
    """test adding k-mers to the genome"""
    gen = def_genome
    gen.add_kmer(25,32)
    assert gen.specific == 25
    assert gen.multimapping == 32
    gen.add_kmer(0,0)
    assert gen.specific == 0
    assert gen.multimapping == 0
    assert gen.length == 9



