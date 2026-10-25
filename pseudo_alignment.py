import pandas as pd
from typing import List,Dict
from read import Read

class Pseudo_alignment:
    def __init__(self, df: pd.DataFrame ,read_number:int ,read_dict: Dict[str,Read], p:int , m:int ):
        """this class is responsible for determining the alignment of the read"""
        self.__df: pd.DataFrame = df
        self.__read_number: int = read_number
        self.__p: int = p
        self.__m: int = m
        self.__skip: List[int] = [0, 0, 0]
        self.__mapping:str = None
        self.__specifique_df: pd.DataFrame = pd.DataFrame(columns=["Genome", "Read", "Read_Position", "Genome_Position"])
        self.__multimapping_df: pd.DataFrame = pd.DataFrame(columns=["Genome", "Read", "Read_Position", "Genome_Position"])
        self.__spe_len: int = 0
        self.__multi_len: int = 0
        self.read_dict: Dict[str, Read] = read_dict
        if not self.df.empty:
            self.__read_apearance = self.df.groupby("Read")["Genome"].nunique()
            self.__specifique_df = df[df["Read"].isin(self.__read_apearance[self.__read_apearance == 1].index)]
            self.__multimapping_df = df[df["Read"].isin(self.__read_apearance[self.__read_apearance > 1].index)]
            self.__spe_len = len(self.__read_apearance[self.__read_apearance == 1].index)
            self.__multi_len = len(self.__read_apearance[self.__read_apearance > 1].index)
            

    @property
    def read_number(self)-> int:
        """this func will return the number of read"""
        return self.__read_number
    
    @property
    def mapping(self)-> str:
        """this func will return if the mapping is unique, ambiguouse ..."""
        return self.__mapping
    
    @property
    def df(self)-> pd.DataFrame:
        """this func will return the panda dataframe"""
        return self.__df
    
    @property
    def specifique_df(self)-> pd.DataFrame:
        """this func will return a panda dataframe with only the specifique reads"""
        return self.__specifique_df
    
    @property
    def multimapping_df(self)-> pd.DataFrame:
        """this func will return a panda dataframe with only the multimapping reads"""
        return self.__multimapping_df
    
    @property
    def spe_len(self)-> int:
        """the nuber of different read in the sepcifiqu df"""
        return self.__spe_len
    
    @property
    def multi_len(self)-> int:
        """the nuber of different read in the multimapping df"""
        return self.__multi_len
    
    @property
    def skip(self)-> List[int]:
        """this func will return the value that were skipped/filtered"""
        return self.__skip
    
    def skip_setter(self,value):
        """Setter for skip"""
        self.__skip = value


    def aligne(self)->List[str]:
        """this function do the alignment algorithme"""
        genome_lst: List[str] = []
        #empty df case or no secific
        if self.df.empty or self.specifique_df.empty:
            self.__mapping = "Unmapped"
            return genome_lst
        
        #if only one specific genome read is uniquely mapped
        unique_genome:str = None
        if len(self.specifique_df["Genome"].unique()) == 1:
            self.__mapping = "uniquely map"
            unique_genome = self.specifique_df["Genome"].unique()[0]
            #one genome mapped case
            if len(self.multimapping_df["Genome"].unique()) == 0:
                return [unique_genome]
             
        #if more than one genome per specific read find the first and the second genome and compare to m
        first_genome = (unique_genome,self.df[self.df["Genome"] == unique_genome]["Read"].count())
        if len(self.specifique_df["Genome"].unique()) > 1:
            total_spe_read_per_genome = self.specifique_df.groupby("Genome")["Read"].nunique().sort_values(ascending=False)
            first_genome = (total_spe_read_per_genome.index[0],total_spe_read_per_genome.iloc[0])#(genome name, number of specifique read)
            second_genome = (total_spe_read_per_genome.index[1],total_spe_read_per_genome.iloc[1])#(genome name, number of specifique read)
            # genome is uniquely map
            if first_genome[1] -second_genome[1] >= self.__m:
                self.__mapping = "uniquely map"
                unique_genome = first_genome[0]
            #first case ambiguously print all genome with specific read
            else:
                self.__mapping = "ambiguously map"
                for genome in self.specifique_df["Genome"].unique():
                    genome_lst.append(genome)
                return genome_lst
        #second part to check by multimapped read
        if self.__mapping == "uniquely map" and self.__p>=0:
            #count the total of read par genome and substract to the previous count of specifique read
            total_read_counts = (self.df.groupby("Genome")["Read"].count().sort_values(ascending=False))
            MaxCount = total_read_counts.iloc[0]
            MapCount =self.df[self.df["Genome"] == first_genome[0]]["Read"].count()

            #second case ambiguously print all of the genome with total count bigger than total count of the first
            if MaxCount - MapCount > self.__p:
                self.__mapping = "ambiguously map"
                kmer_per_genome= self.df.groupby("Genome")["Read"].count()
                for genome,count in kmer_per_genome[kmer_per_genome >= MapCount].items():
                    genome_lst.append(genome)
                return genome_lst
            else:
                return [unique_genome]

    def filter_by_read_apearance(self,apearance:int)->int:
        """this function filter the df by the number of genome a kmer apear and return the skipped number"""
        genome_counts = self.df.groupby(["Read", "Read_Position"])["Genome"].nunique()
        reads_to_remove = genome_counts[genome_counts > apearance].index
        filtered_df = self.df[~self.df.set_index(["Read", "Read_Position"]).index.isin(reads_to_remove)]
        count = len(self.df)- len(filtered_df)
        self.__read_number -= count
        self.__df = filtered_df
        self.__read_apearance = self.df.groupby("Read")["Genome"].nunique()
        self.__specifique_df = self.df[self.df["Read"].isin(self.__read_apearance[self.__read_apearance == 1].index)]
        self.__multimapping_df = self.df[self.df["Read"].isin(self.__read_apearance[self.__read_apearance > 1].index)]
        self.__spe_len = len(self.__read_apearance[self.__read_apearance == 1].index)
        self.__multi_len = len(self.__read_apearance[self.__read_apearance > 1].index)
        return count
    
    def filter_by_read_quality_of_kmer(self,quality:str,k:int)->int:
        newdf = pd.DataFrame(columns=["Genome", "Read", "Read_Position", "Genome_Position"])
        count = 0
        for col, row in self.df.iterrows():
            read_name = row["Read"]
            position = row["Read_Position"]
            quality_seq = self.read_dict[read_name].quality[position:position+k]
            if self.read_dict[read_name].quality_for_kmer(quality_seq) > quality:
                newdf = pd.concat([newdf, pd.DataFrame([row])], ignore_index=True)
            else:
                count+=1
        self.__df = newdf
        self.__read_apearance = self.df.groupby("Read")["Genome"].nunique()
        self.__specifique_df = self.df[self.df["Read"].isin(self.__read_apearance[self.__read_apearance == 1].index)]
        self.__multimapping_df = self.df[self.df["Read"].isin(self.__read_apearance[self.__read_apearance > 1].index)]
        self.__spe_len = len(self.__read_apearance[self.__read_apearance == 1].index)
        self.__multi_len = len(self.__read_apearance[self.__read_apearance > 1].index)
        return count
            