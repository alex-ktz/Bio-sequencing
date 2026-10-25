from typing import Dict,Tuple,List
from read import Read
import numpy as np
import io

class Read_reference:
    def __init__(self, fq_file: str,k: int):
        """this func create a reade reference"""
        self.__fq_file: str = fq_file
        self.__k: int = k
        self.__read_dictionary: Dict[str, Read] = self.read_extractor()
        self.__read_number: int = len(self.read_dictionary)
        self.__read_array: Tuple[List[np.ndarray], List[int]] = self.reads_dict_to_array()
        

    @property
    def read_number(self)-> int:
        """this func will return the number of read"""
        return self.__read_number
      
    @property
    def read_dictionary(self)-> Dict[str,Read]:
        """this func will return the read dictionary"""
        return self.__read_dictionary
    
    @property
    def read_array(self)-> np.ndarray:
        """this func will return the read array"""
        return self.__read_array
    
    def __str__(self)-> str:
        """this func will print the read ref"""
        return f"there are {self.read_number} reads in the refference"
    
    def read_extractor(self)-> Dict[str, Read]:
        """this function recieve a fq file and return a dict of reads {read id: Read object}"""
        read_dict: Dict[str,Read] = {}
        try:
            with io.open(self.__fq_file, 'r') as file:
                while True:
                    #read the next 4 lines at once
                    lines = [file.readline().strip() for _ in range(4)]
                    #if the first line is empty, we've reached the end of the file
                    if not lines[0]:
                        break
                    at = lines[0][:1]
                    read_id = lines[0][1:]  # skip @
                    read_sequence = lines[1]
                    plus = lines[2]
                    read_quality = lines[3]
                    #if one line is missing
                    if at!= "@" or plus !="+":
                        print("error in a read")
                        if read_quality[:1] == "@":
                            errorlines = [file.readline().strip() for _ in range(3)]
                            read_id = read_quality[1:]
                            read_sequence = errorlines[0]
                            plus = errorlines[1]
                            read_quality = errorlines[2]
                    # Validate the read
                    if read_sequence and len(read_sequence) == len(read_quality) and read_id and set(read_sequence).issubset({"A", "C", "G", "T", "N"}) and plus =="+" and at == "@":
                        #check for doubles
                        read_dict.setdefault(read_id, Read(read_sequence, read_id, read_quality))
                    else:
                        print("error in a read")
        except FileNotFoundError:
            print("Error: The FASTQ file was not found.")
        except UnicodeDecodeError:
            print("Error: The file contains invalid characters or is not properly encoded.")
        except ValueError as e:
            print(f"Error in FASTQ file: {e}")
        except Exception as e:
            print(f"General error in FASTQ file: {e}")
        return read_dict
    
    @staticmethod
    def convert_string_to_np(sequence: str)-> np.ndarray:
        """convert an sequence of A,C,G,T,N to A=0, C=1, G=2, T=3 N=4 and return it as a numpy array to improve the efficiency of searching"""
        mapping = np.array([-1] * 256, dtype=np.int8)
        mapping[ord('A')] = 0
        mapping[ord('C')] = 1
        mapping[ord('G')] = 2
        mapping[ord('T')] = 3
        mapping[ord('N')] = 4
        # convert sequence to a numpy array of ASCII values of the letter i  one vectorised step
        ascii_arr:np.ndarray = np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)
        # using the mapping built earlier and apply it on the ascii array (wich will happen in O(1))
        mapped_seq:np.ndarray = mapping[ascii_arr]
        # check for invalid characters
        if np.any(mapped_seq == -1):
            raise ValueError("invalid character in sequence")
        return mapped_seq

    def reads_dict_to_array(self) -> Tuple[List[np.ndarray], List[int]]:
        """convert all the reads into k-mer arrays of number(with a list that keep the read index)"""
        read_kmers: List[np.ndarray] = []
        lst_of_index: List[int] = []
        for read_idx, (read_name, read_obj) in enumerate(self.__read_dictionary.items()):
            read_seq = Read_reference.convert_string_to_np(read_obj.sequence)
            if len(read_seq) < self.__k:
                continue  # Skip reads shorter than k
            read_kmer_array = np.lib.stride_tricks.sliding_window_view(read_seq, self.__k)
            read_kmers.append(read_kmer_array)
            # Store read index to retrieve later
            lst_of_index.append(read_idx)
        return read_kmers, lst_of_index
    

    def reads_dict_to_array_complement(self)-> np.ndarray:
        """convert all the *complment* reads into k-mer arrays of number(with a list that keep the read index)"""
        read_kmers: List[np.ndarray] = []
        lst_of_index: List[int] = []
        for read_idx, (read_name, read_obj) in enumerate(self.__read_dictionary.items()):
            read_seq = Read_reference.convert_string_to_np(read_obj.complement)
            if len(read_seq) < self.__k:
                continue  # Skip reads shorter than k
            read_kmer_array = np.lib.stride_tricks.sliding_window_view(read_seq, self.__k)
            read_kmers.append(read_kmer_array)
            # Store read index to retrieve later
            lst_of_index.append(read_idx)
        return read_kmers, lst_of_index

    def filter_by_read_quality(self,filter_value)-> int:
        """this func delete all of the kmer that were not of a good enough quality while keeping a track of how many it deleted"""
        new_dict: Dict[str, Read] = {}
        count = 0
        for read_id, read_obj in self.read_dictionary.items():
            if read_obj.quality_score > filter_value:
                new_dict[read_id] = read_obj
            else:
                count += 1
        self.__read_dictionary = new_dict
        self.__read_array = self.reads_dict_to_array()
        self.__read_number = len(self.read_dictionary)
        return count
    
