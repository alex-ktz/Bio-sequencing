import pytest
from read_reference import Read_reference
import numpy as np

@pytest.fixture
def sample_read_ref():
    """create default sample k-mers, genomes, and numpy arrays for testing."""
    broken_read = Read_reference("test\\broken kmer.fq",4)
    ten_read = Read_reference("test\\ten kmer.fq",4)
    two_read =Read_reference("test\\two kmer.fq",2)
    emmpty_read =Read_reference("test\\empty kmer.fq",4)

    return broken_read, ten_read, two_read, emmpty_read

### Tests ###
def test_read_extractor(sample_read_ref):
    """test that read properly creates a genome dictionary with the value in the FASTQ even with complicated value"""
    broken_read, ten_read, two_read, emmpty_read= sample_read_ref
    assert len(ten_read.read_dictionary) == 10
    assert ten_read.read_number == 10
    assert "kmer1" in ten_read.read_dictionary
    assert "kmer2" in ten_read.read_dictionary
    assert "kmer3" in ten_read.read_dictionary
    assert ten_read.read_dictionary["kmer3"].sequence == "GGGCGTCCAA"
    assert ten_read.read_dictionary["kmer7"].sequence == "GAACAAGAGC"
    assert ten_read.read_dictionary["kmer8"].sequence == "AAGACGCCAA"
    assert ten_read.read_dictionary["kmer5"].quality == "?DE?EBB@CD"
    assert ten_read.read_dictionary["kmer2"].quality == "DBCC?>??@@"
    assert "correct_kmer" in broken_read.read_dictionary
    assert broken_read.read_number == 1 
    assert len(emmpty_read.read_dictionary) == 0
    assert emmpty_read.read_number == 0

    
def test_convert_string_to_np():
    """Test conversion of a DNA string to a numpy array."""
    seq = "ACGTN"
    expected = np.array([0, 1, 2, 3, 4], dtype=np.int8)
    result = Read_reference.convert_string_to_np(seq)
    np.testing.assert_array_equal(result, expected)
    seq = ""
    expected = np.array([], dtype=np.int8)
    result = Read_reference.convert_string_to_np(seq)
    np.testing.assert_array_equal(result, expected)
    





