from typing import Dict,List,Optional,Any
from genome import Genome
import numpy as np
import io

class Genome_reference:
    def __init__(self, fa_file: str, k: int):
        """this func create a genome refference with a genome dict to store all the data from each genome and a genome ref to dump"""
        self.__fa_file: str = fa_file
        self.__k: int = int(k)
        self.__genome_dictionary: Dict[str, Genome] = self.genome_extractor()
        self.__kmer_ref: Optional[Dict] = None
        self.__similarity_dict: Optional[Dict] = None

    @property
    def genome_dictionary(self)-> Dict[str,Genome]:
        """return the genome dict"""
        return self.__genome_dictionary

    @property
    def k(self)-> List[np.ndarray]:
        """return the k"""
        return self.__k
    
    @property
    def similarity_dict(self)-> Dict[str,Genome]:
        """return the genome dict"""
        return self.__similarity_dict
    
    def __str__(self)-> str:
        """this func will print the genome ref"""
        return f"there are {len(self.__genome_dictionary)} genomes in the refference"
    
    def genome_extractor(self)-> Dict[str, Genome]:
            """this function recieve a fa file and return a dict of genome {genome id: Genome object}"""
            #chek for unvalid char
            genome_dict:Dict[str,Genome] = {}
            genome_id: str = None
            genome_sequence: List[str] = []
            valid_nucleotides = {"A", "C", "G", "T", "N"}
            try:
                with io.open(self.__fa_file, 'r', encoding="utf-8") as file:
                    for line in file:
                        line = line.strip()
                        if line.startswith(">"):
                            if genome_id is not None:# if the genome is not the first one enter the last Genome in the dict
                                seq = "".join(genome_sequence)
                                if seq and set(seq).issubset(valid_nucleotides):
                                    genome_dict[genome_id] = Genome(genome_id, seq)
                                else:
                                    print("error in a genome")
                            genome_id = line[1:]  # skip the >
                            genome_sequence = []# empty the seq to recieve the next one
                        else:
                            genome_sequence.append(line)# add that line to the sequence
                    #save the last genome
                    if genome_id:
                        if len("".join(genome_sequence))>0 and set("".join(genome_sequence)).issubset({"A", "C", "G", "T", "N"}):#checking validity of the character in O(1)
                            genome_dict[genome_id] = Genome(genome_id, "".join(genome_sequence))
                        seq = "".join(genome_sequence)
                        if seq and set(seq).issubset(valid_nucleotides):
                            genome_dict[genome_id] = Genome(genome_id, seq)
                        else:
                            print("Error in genome")
            except FileNotFoundError:
                print("Error: The file was not found.")
            except UnicodeDecodeError:
                print("Error: The file contains invalid characters or is not properly encoded.")
            except ValueError as e:
                print(f"Error in FASTA file: {e}")
            except Exception as e:
                print(f"General error in FASTA file: {e}")
            return genome_dict

    @staticmethod
    def convert_string_to_np(sequence: str)-> np.ndarray:
        """convert an sequence of A,C,G,T,N to A=0, C=1, G=2, T=3 N=4 and return it as a numpy array to improve the efficiency of searching"""
        mapping :np.ndarray= np.array([-1] * 256, dtype=np.int8)
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

    def genome_dict_to_lst_of_array(self) -> List[np.ndarray]:
        """the function recieve a list of genome array and return an nparray of all of the possible kmer in the genome"""
        genome_lst: List[np.ndarray] = []
        for id,genome in self.__genome_dictionary.items():
            genome_arr = Genome_reference.convert_string_to_np(genome.sequence)
            if self.__k > genome.length:
                raise ValueError("k must be smaller than the genome")
            # numpyh stride is used for fast windowing and is a lot more efficient than doing a manual loop
            genome_lst.append(np.lib.stride_tricks.sliding_window_view(genome_arr, window_shape = self.__k))
        return genome_lst

    def create_kmer_ref(self)-> Dict[str,Dict[str, List[int]]]:
        """this func create the dict that is used to dump the ref of the genome"""
        kmer_occurrences: Dict[str,Dict[str, List[int]]] = {}
        for genome_id, genome_obj in self.genome_dictionary.items():
            seq: str= genome_obj.sequence
            seq_len: int = len(seq)
            for pos in range(seq_len - self.__k + 1):  # Sliding window
                kmer:str = seq[pos : pos + self.__k]
                genome_positions = kmer_occurrences.setdefault(kmer, {})
                genome_positions.setdefault(genome_id, []).append(pos)
        self.__kmer_ref = kmer_occurrences
        return kmer_occurrences
    
    def kmer_per_genome(self)-> None:
        """this func add to each Genome object of the refference a value for their unique and their multimapping kmer"""
        gendict: Dict[str,List[int]] = {}
        #checkin if kmer ref already exist to avoid double work
        kmer_ref = self.__kmer_ref if self.__kmer_ref is not None else self.create_kmer_ref()
        for kmer, genomes in kmer_ref.items():
            is_specific:bool = len(genomes) > 1 
            for genome_name in genomes:
                if genome_name not in gendict:
                    gendict[genome_name] = [0, 0] 
                if is_specific:
                    gendict[genome_name][1] += 1  
                else:
                    gendict[genome_name][0] += 1 
        for genome in self.genome_dictionary.items():
            spe_kmer = gendict[genome[0]][0]
            multimapping_kmer = gendict[genome[0]][1]
            genome[1].add_kmer(spe_kmer,multimapping_kmer)
  
    def sorting_for_similarity(self)-> None:
        """sort the genomes by the requested order"""
        #making sure the genome have thei kmer
        self.kmer_per_genome()
        kmer_ref = self.__kmer_ref if self.__kmer_ref is not None else self.create_kmer_ref()
        genome_list = list(self.genome_dictionary.values())
        sorted_genomes:List[Genome] = sorted(genome_list, key=lambda g: (g.specific, g.multimapping, g.length), reverse=True)
        return sorted_genomes

    def filter_by_similarity(self,filter_value:int):
        """filter the genomes by similarity"""
        sorted_list: List[Genome] = self.sorting_for_similarity()
        similarity_dict: Dict[str, Dict[str, Any]] = {}
        #building the dict that will be dumped
        for genome in sorted_list:
            similarity_dict[genome.id] = {"kept":"yes", "unique_kmers":genome.specific , "total_kmers":genome.specific+genome.multimapping, "genome_length":genome.length, "similar_to":"NA", "similarity_score":"NA"}
        #putting all of the kmer of both genome in a set
        for i in range(len(sorted_list)):
            for j in range(i + 1, len(sorted_list)):
                kmer_set = set()
                for kmer,genome_list in self.__kmer_ref.items():
                        if sorted_list[i].id in genome_list or sorted_list[j].id in genome_list:
                            kmer_set.add(kmer)
                if kmer_set:
                    minimum = min(sorted_list[i].length,sorted_list[j].length)
                    similarity_ratio = minimum / len(kmer_set)
                    #updating the dict
                    if similarity_ratio > filter_value:
                        similarity_dict[sorted_list[j].id]["kept"] ="no"
                        similarity_dict[sorted_list[j].id]["similar_to"] =sorted_list[i].id
                        similarity_dict[sorted_list[j].id]["similarity_score"] = (minimum/len(kmer_set))
        self.__similarity_dict = similarity_dict
        #deleting in the dict
        for genome_id, info in list(similarity_dict.items()):
            if info["kept"] == "no":
                del self.__genome_dictionary[genome_id]
