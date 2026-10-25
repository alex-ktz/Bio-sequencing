from pseudo_alignment import Pseudo_alignment
import pandas as pd

def test_align():
    """testing the align algo"""
    #empty df
    data = []
    dataf = Pseudo_alignment(pd.DataFrame(data),0,{},1,1)
    assert dataf.aligne() == []

    #df with only one genome
    data = [{"Genome": "a", "Read": "kmer1", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),1,{},1,1)
    assert dataf.aligne() == ['a']

    # two genome with different kmer
    data = [{"Genome": "a", "Read": "kmer1", "Genome_position": 176, "Read_position":48},{"Genome": "b", "Read": "kmer2", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),2,{},1,1)
    assert dataf.aligne() == ['a','b']

    # two genome with same kmer (unmapped because no specific)
    data = [{"Genome": "a", "Read": "kmer1", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer1", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),2,{},1,1)
    assert dataf.aligne() == []

    # two genome with a difference in specifique bigger than m
    data = [{"Genome": "a", "Read": "kmer4", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer1", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer2", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer3", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),4,{},1,1)
    assert dataf.aligne() == ['b']

    # two genome with a difference in specifique smaller than m ambi to all spe
    data = [{"Genome": "a", "Read": "kmer1", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer2", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer3", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer4", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),4,{},1,4)
    assert dataf.aligne() == ['a','b']

    # three genome with a difference in specifique smaller than m ambi to all spe
    data = [{"Genome": "a", "Read": "kmer1", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer2", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer3", "Genome_position": 176," Read_position":48},{"Genome": "c", "Read": "kmer4", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),4,{},1,4)
    assert dataf.aligne() == ['a','b','c']

    # three genome with a difference in specifique bigger than m
    data = [{"Genome": "a", "KmReader": "kmer1", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer2", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer3", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer4", "Genome_position": 176," Read_position":48},{"Genome": "c", "Read": "kmer5", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),5,{},1,1)
    assert dataf.aligne() == ['b']

    # # three genome with a difference in specifique smaller than m ambi to all spe but only two spe
    data = [{"Genome": "a", "Read": "kmer1", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer2", "Genome_position": 176," Read_position":48},{"Genome": "b", "Read": "kmer3", "Genome_position": 176," Read_position":48},{"Genome": "c", "Read": "kmer3", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),4,{},1,4)
    assert dataf.aligne() == ['a','b']

# max - map > P so ambi read with all the bigger than map (b)
    data = [
        {"Genome": "a", "Read": "kmer1", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer2", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer3", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer4", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer5", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer6", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer7", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer5", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer8", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer3", "Genome_position": 176," Read_position":48},
        {"Genome": "c", "Read": "kmer1", "Genome_position": 176," Read_position":48},
        {"Genome": "c", "Read": "kmer2", "Genome_position": 176," Read_position":48},
        {"Genome": "c", "Read": "kmer4", "Genome_position": 176," Read_position":48}]
    dataf = Pseudo_alignment(pd.DataFrame(data),1,{},1,1)
    assert dataf.aligne() == ['a','b']

    #max - map < P so unique read
    data = [
        {"Genome": "a", "Read": "kmer1", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer2", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer3", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer4", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer5", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer6", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer7", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer8", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer9", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer4", "Genome_position": 176," Read_position":48},
        {"Genome": "c", "Read": "kmer1", "Genome_position": 176," Read_position":48},
        {"Genome": "c", "Read": "kmer3", "Genome_position": 176," Read_position":48},
        {"Genome": "c", "Read": "kmer5", "Genome_position": 176," Read_position":48},]
    dataf = Pseudo_alignment(pd.DataFrame(data),3,{},3,1)
    assert dataf.aligne() == ['b']

    # ambi read to all spe
    data = [
        {"Genome": "a", "Read": "kmer1", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer2", "Genome_position": 176," Read_position":48},
        {"Genome": "a", "Read": "kmer3", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer1", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer4", "Genome_position": 176," Read_position":48},
        {"Genome": "b", "Read": "kmer5", "Genome_position": 176," Read_position":48},
        {"Genome": "c", "Read": "kmer6", "Genome_position": 176," Read_position":48},
        {"Genome": "c", "Read": "kmer7", "Genome_position": 176," Read_position":48},
        {"Genome": "d", "Read": "kmer5", "Genome_position": 176," Read_position":48},
        {"Genome": "f", "Read": "kmer1", "Genome_position": 176," Read_position":48}
        ]
    dataf = Pseudo_alignment(pd.DataFrame(data),10,{},1,1)
    assert dataf.aligne() == ['a', 'b', 'c']
        