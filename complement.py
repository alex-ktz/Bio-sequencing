from pseudo_alignment import Pseudo_alignment


class Complement:
    def __init__(self,normal_align: Pseudo_alignment,complement_align: Pseudo_alignment):
        self.__normal_align: Pseudo_alignment = normal_align
        self.__complement_align: Pseudo_alignment = complement_align
        self.__winner: Pseudo_alignment = None
        self.__looser: Pseudo_alignment = None
        self.choose_read()
    
    @property
    def winner(self)-> Pseudo_alignment:
        """this func will return only the winner align"""
        return self.__winner
    
    @property
    def looser(self)-> Pseudo_alignment:
        """this func will return only the looser align"""
        return self.__looser

    def choose_read(self)->None:
        self.__normal_align.aligne()
        self.__complement_align.aligne()
        normal_mapping: str = self.__normal_align.mapping
        comp_mapping: str = self.__complement_align.mapping
        if normal_mapping == "uniquely map" and comp_mapping != "uniquely map":
            self.__winner = self.__normal_align
            self.__looser = self.__complement_align
        elif normal_mapping != "uniquely map" and comp_mapping == "uniquely map":
            self.__winner = self.__complement_align
            self.__looser = self.__normal_align

        elif normal_mapping == "uniquely map" and comp_mapping == "uniquely map":
            if len(self.__normal_align.specifique_df) >= len(self.__complement_align.specifique_df):
                self.__winner = self.__normal_align
                self.__looser = self.__complement_align
            else:
                self.__winner = self.__complement_align
                self.__looser = self.__normal_align

        elif normal_mapping == "ambiguously map" and comp_mapping == "ambiguously map":
            if len(self.__normal_align.df) >= len(self.__complement_align.df):
                self.__winner = self.__normal_align
                self.__looser = self.__complement_align
            else:
                self.__winner = self.__complement_align
                self.__looser = self.__normal_align
        else:
            self.__winner = self.__normal_align
            self.__looser = self.__complement_align
        

