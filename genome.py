from typing import Optional
class Genome:
    def __init__(self, id: str, sequence: str) -> None:
        """Create a new Genome object."""
        self.__id: str = id
        self.__sequence: str = sequence
        self.__length: int = len(sequence)
        self.__specific: Optional[int] = None
        self.__multimapping: Optional[int] = None

    @property
    def id(self)-> str:
        """this func will return the id of the genome"""
        return self.__id
    @property
    def sequence(self)-> str:
        """this func will return the full sequence of the genome"""
        return self.__sequence
    
    @property
    def length(self)-> int:
        """this func will return the length of the genome"""
        return self.__length
    
    @property
    def specific(self)-> int:
        """this func will return the nuber of specific kmer of the genome"""
        return self.__specific
    
    @property
    def multimapping(self)-> int:
        """this func will return the multimapping kmer of the genome"""
        return self.__multimapping
    
    def __str__(self)-> str:
        """this func will print the genome"""
        return f"id of the genome: {self.__id}   specific kmer {self.__specific}  multimapping kmer {self.__multimapping}" 

    def add_kmer(self, spe: int, multimapping: int) -> None:
        """this func add to the genome his number of specific and multimapping kmer"""
        self.__specific = spe
        self.__multimapping = multimapping