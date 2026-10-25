from read import Read
import pytest

@pytest.fixture
def def_read():
    """Making a default read for the basic case"""
    return Read("ATCG","girafe","QUAL")

def test_opposite(def_read):
    """Test the get opposite func"""
    assert def_read.complement == "CGAT"
    #invalid value
    assert def_read.complement != "GCTA"
    #trying the N value
    read = Read("ATNCG","girafe","ABCDE")
    assert read.complement == "CGNAT"

def test_quality_for_kmer(def_read)-> float:
        """test quality for kmer"""
        assert def_read.quality_score == 43.75
        assert def_read.quality_for_kmer("cdetgr)") == 63.285714285714285



