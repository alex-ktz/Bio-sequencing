from typing import Dict,List

class Read:
    def __init__(self,sequence: str, id: str, quality: str):
        """this func create a new read with its sequence, id and quality string"""
        self.__sequence: str = sequence
        self.__id: str = id
        self.__quality: str = quality
        

    @property
    def sequence(self)-> str:
        """this func will return the sequence of the read as a string"""
        return self.__sequence
    
    @property
    def id(self) -> str:
        """this func will return the id of the read"""
        return self.__id
    
    @property
    def quality(self)-> str:
        """this func will return the quality of the sequence of the read"""
        return self.__quality
    
    @property
    def quality_score(self)-> float:
        """this func will return the quality score of that read"""
        return self.quality_for_kmer(self.__quality)
    
    @property
    def complement(self)-> str:
        """this func will return the opposite strand of the read"""
        return self.complement_DNA()  

    def complement_DNA(self)-> str:
        """this func give the reverse complement strand of the DNA like it will appear in the genome"""
        complement_map: Dict[str, str] = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        if not set(self.__sequence).issubset({"A", "C", "G", "T", "N"}):
            print("Invalid letter in read sequence")
        return "".join(complement_map[nuc] for nuc in self.__sequence)[::-1]
    
    def __str__(self)-> str:
        """this func will print the read"""
        return f"id of the read: {self.__id}\nsequence: {self.__sequence}\nopposite: {self.complement}\nquality : { self.__quality}"

    def quality_for_kmer(self,seq)-> float:
        """this func will return the quality score of a given kmer"""
        phred_scores: List[int] = [ord(q) - 33 for q in seq]
        return sum(phred_scores) / len(phred_scores) 


