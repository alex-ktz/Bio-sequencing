from typing import Dict,List,Tuple
import numpy as np
import pandas as pd
from read import Read
from genome import Genome


class Matching_algorythme:
    def __init__(self, read_array: Tuple[List[np.ndarray], List[int]], read_dictionary: Dict[str,Read], lst_of_genome_array: List[np.ndarray], genome_dictionary: Dict[str,Genome]):
        """this func create a matching algorythme that will give us the panda dataframe well need for the alignment"""
        self.__read_lst: List[str] = list(read_dictionary.keys())
        self.__genome_lst: List[str] = list(genome_dictionary.keys())
        self.__read_array: Tuple[List[np.ndarray], List[int]] = read_array
        self.__lst_of_genome_array: List[np.ndarray] = lst_of_genome_array
        self.__matches_array: np.ndarray = self.match_reads_in_genome()
        self.__panda_data_frame: pd.DataFrame = self.matches_to_dataframe()
    
    @property
    def panda_data_frame(self)-> pd.DataFrame:
        """return tha panda dataframe"""
        return self.__panda_data_frame


    def match_reads_in_genome(self) -> np.ndarray:
        """Match k-mers from reads against k-mers in genome and return matches [genome index, read index, read k-mer position, genome position]."""
        read_kmers, read_indexes = self.__read_array
        # Dictionary for fast lookup {kmer_tuple: [(read_idx, read_pos)]}
        read_dict: Dict[Tuple[int], List[Tuple[int, int]]] = {}
        for read_idx, read_kmer_array in zip(read_indexes, read_kmers):
            for read_pos, kmer in enumerate(read_kmer_array):
                kmer_tuple = tuple(kmer)
                #create a dict of different kmer with thei index and position
                if kmer_tuple not in read_dict:
                    read_dict[kmer_tuple] = []
                read_dict[kmer_tuple].append((read_idx, read_pos))  # store read index and k-mer position
        matches = []
        for genome_idx, genome_array in enumerate(self.__lst_of_genome_array):
            for genome_pos, genome_kmer in enumerate(genome_array):
                kmer_tuple = tuple(genome_kmer)
                if kmer_tuple in read_dict:
                    for read_idx, read_pos in read_dict[kmer_tuple]:
                        matches.append((genome_idx, read_idx, read_pos, genome_pos))
        return np.array(matches, dtype=np.int32)


    def matches_to_dataframe(self) -> pd.DataFrame:
        """transform the np.array of matches into a pandas DataFrame for later analysis."""
        records = []
        for genome_idx, read_idx, read_pos, genome_pos in self.__matches_array:
            data = {"Genome": self.__genome_lst[genome_idx], "Read": self.__read_lst[read_idx], "Read_Position": read_pos,"Genome_Position": genome_pos}
            records.append(data)
        return pd.DataFrame(records)














