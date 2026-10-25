import pytest
import numpy as np
import pandas as pd
from matching_algorythme import Matching_algorythme
from read import Kmer
from genome import Genome

@pytest.fixture
def sample_data():
    """create default sample k-mers, genomes, and numpy arrays for testing."""
    kmer_dict = {"kmer1": Kmer("AAA", "kmer1", "###"),"kmer2": Kmer("TTT", "kmer2", "###"),"kmer3": Kmer("CCC", "kmer3", "###")}
    genome_dict = {"Genome_1": Genome("Genome_1", "AAATTTCCC"),"Genome_2": Genome("Genome_2", "TTTCCCAAA")}
    kmer_array = np.array([[0, 0, 0],[3, 3, 3],[1, 1, 1]], dtype=np.int8)
    genome_array_list = [np.array([[0, 0, 0],[0, 0, 3],[0, 3, 3],[3, 3, 3],[3, 3, 1],[3, 1, 1],[1, 1, 1]], dtype=np.int8), np.array([[3, 3, 3],[3, 3, 1],[3, 1, 1],[1, 1, 1],[1, 1, 0],[1, 0, 0],[0, 0, 0]], dtype=np.int8)]

    return kmer_array, kmer_dict, genome_array_list, genome_dict


def test_match_kmers_in_genome(sample_data):
    """test if k-mers are correctly matched in the genome"""
    kmer_array, kmer_dict, genome_array_list, genome_dict = sample_data
    match_algo = Matching_algorythme(kmer_array, kmer_dict, genome_array_list, genome_dict)
    matches = match_algo.match_kmers_in_genome()
    expected_matches = np.array([[0, 0, 0],[0, 1, 3],[0, 2, 6],[1, 1, 0],[1, 2, 3],[1, 0, 6]], dtype=np.int32)
    np.testing.assert_array_equal(matches, expected_matches)
    #empty genome array
    match_algo = Matching_algorythme(kmer_array, kmer_dict, np.array(([])), genome_dict)
    matches = match_algo.match_kmers_in_genome()
    expected_matches = np.array(([]))
    np.testing.assert_array_equal(matches, expected_matches)
    #empty kmer array
    match_algo = Matching_algorythme(np.array(([])), kmer_dict, genome_array_list, genome_dict)
    matches = match_algo.match_kmers_in_genome()
    expected_matches = np.array(([]))
    np.testing.assert_array_equal(matches, expected_matches)



def test_matches_to_dataframe(sample_data):
    """check if the matches are correctly converted to a dataframe"""
    kmer_array, kmer_dict, genome_array_list, genome_dict = sample_data
    match_algo = Matching_algorythme(kmer_array, kmer_dict, genome_array_list, genome_dict)
    actual_df = match_algo.matches_to_dataframe()
    # make sure the DataFrame has the correct columns
    assert list(actual_df.columns) == ["Genome", "Kmer", "Position"]
    # check if the df has expected values
    expected_data = {"Genome": ["Genome_1", "Genome_1", "Genome_1", "Genome_2", "Genome_2", "Genome_2"],"Kmer":   ["kmer1",     "kmer2",     "kmer3",     "kmer2",     "kmer3",     "kmer1"],"Position": [0, 3, 6, 0, 3, 6]}
    expected_df = pd.DataFrame(expected_data)
    pd.testing.assert_frame_equal(actual_df, expected_df,check_dtype=False)
    #empty set
    match_algo = Matching_algorythme(np.array(([])), kmer_dict, genome_array_list, genome_dict)
    actual_df = match_algo.matches_to_dataframe()
    # make sure the DataFrame has the correct columns
    assert list(actual_df.columns) == []
    # check if the df has expected values
    expected_data = {}
    expected_df = pd.DataFrame(expected_data)
    pd.testing.assert_frame_equal(actual_df, expected_df,check_dtype=False)
