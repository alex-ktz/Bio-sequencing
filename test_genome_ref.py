import pytest
import numpy as np
from genome_reference import Genome_reference

@pytest.fixture
def sample_genome_ref():
    """create default sample k-mers, genomes, and numpy arrays for testing."""
    one_gen:Genome_reference = Genome_reference("test\\one genome.fa", 3)
    three_gen:Genome_reference = Genome_reference("test\\.three genome.fa", 3)
    broken_gen:Genome_reference =Genome_reference("test\\uncorrect genome.fa", 3)
    empty_gen:Genome_reference = Genome_reference("test\\empty genome.fa", 3)
    two_gen:Genome_reference = Genome_reference("test\\two genome.fa", 3)

    return three_gen, one_gen, broken_gen, empty_gen, two_gen

### Tests ###
def test_genome_extractor(sample_genome_ref):
    """test that genome_extractor properly creates a genome dictionary with the value in the FASTA"""
    three_gen, one_gen, broken_gen, empty_gen, two_gen = sample_genome_ref
    assert len(three_gen.genome_dictionary) == 3
    assert "genome1" in three_gen.genome_dictionary
    assert "genome2" in three_gen.genome_dictionary
    assert "genome3" in three_gen.genome_dictionary
    assert three_gen.genome_dictionary["genome1"].sequence == "AAACCCGGGTTT"
    assert three_gen.genome_dictionary["genome2"].sequence == "TTTAAATTTGGG"
    assert three_gen.genome_dictionary["genome3"].sequence == "TTTAAATTTGGA"
    assert len(one_gen.genome_dictionary) == 1
    assert "lonely_genome" in one_gen.genome_dictionary
    assert one_gen.genome_dictionary["lonely_genome"].sequence == "AAACCCGGGTTT"
    assert len(empty_gen.genome_dictionary) == 0
    assert len(broken_gen.genome_dictionary) == 1
    assert broken_gen.genome_dictionary["valid_genome"].sequence == "CTCGATGGGTTTGAGGCGGCGCTGACCGACAAAGCGGCCAAGCAGCGCATGTGGGCCTAC"

def test_convert_string_to_np():
    """Test conversion of a DNA string to a numpy array."""
    seq = "ACGTN"
    expected = np.array([0, 1, 2, 3, 4], dtype=np.int8)
    result = Genome_reference.convert_string_to_np(seq)
    np.testing.assert_array_equal(result, expected)
    seq = ""
    expected = np.array([], dtype=np.int8)
    result = Genome_reference.convert_string_to_np(seq)
    np.testing.assert_array_equal(result, expected)

def test_genome_dict_to_lst_of_array(sample_genome_ref):
    """Test that genome_dict_to_lst_of_array returns the correct data"""
    three_gen, one_gen, broken_gen, empty_gen, two_gen = sample_genome_ref
    arr = one_gen.genome_dict_to_lst_of_array()
    expected = np.array(([[[0, 0, 0],[0, 0, 1],[0, 1, 1],[1, 1, 1],[1, 1, 2],[1, 2, 2],[2, 2, 2],[2, 2, 3],[2, 3, 3],[3, 3, 3]]]),dtype=np.int8)
    np.testing.assert_array_equal(arr, expected)
    arr = empty_gen.genome_dict_to_lst_of_array()
    expected = np.array(([]))
    np.testing.assert_array_equal(arr, expected)
    arr = two_gen.genome_dict_to_lst_of_array()
    expected = np.array([[[0, 0, 1], [0, 1, 3], [1, 3, 0]], [[3, 3, 2], [3, 2, 1], [2, 1, 0]]], dtype=np.int8)
    np.testing.assert_array_equal(arr, expected)

def test_create_kmer_ref(sample_genome_ref):
    """Test that create_kmer_ref creates a valid k-mer reference dictionary."""
    three_gen, one_gen, broken_gen, empty_gen, two_gen = sample_genome_ref
    assert two_gen.create_kmer_ref() == {'AAC': {'genome1': [0]}, 'ACT': {'genome1': [1]}, 'CTA': {'genome1': [2]}, 'TTG': {'genome2': [0]}, 'TGC': {'genome2': [1]}, 'GCA': {'genome2': [2]}}
    assert empty_gen.create_kmer_ref() =={}

def test_kmer_per_genome(sample_genome_ref):
    """Test that kmer_per_genome properly updates each Genome object with k-mer counts."""
    three_gen, one_gen, broken_gen, empty_gen, two_gen = sample_genome_ref
    three_gen.kmer_per_genome()
    assert three_gen.genome_dictionary['genome1'].specific == 7
    assert three_gen.genome_dictionary['genome1'].multimapping == 3
    assert three_gen.genome_dictionary['genome2'].specific == 0
    assert three_gen.genome_dictionary['genome2'].multimapping == 9
    assert three_gen.genome_dictionary['genome3'].specific == 1
    assert three_gen.genome_dictionary['genome3'].multimapping == 8

def test_filter_by_similarity(sample_genome_ref):
    three_gen, one_gen, broken_gen, empty_gen, two_gen = sample_genome_ref
    three_gen.filter_by_similarity(1.1)
    assert len(three_gen.genome_dictionary) == 2
    
    

