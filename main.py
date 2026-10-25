###########################################################
#alexander katz 332391077
#user-alexkatz 
#mail-alexandr.katz@mail.huji.ac.il
#aditionale notes-None
###########################################################
import argparse
import json
import pickle
from typing import Optional, Any
import pandas as pd
import time
from matching_algorythme import Matching_algorythme
from pseudo_alignment import Pseudo_alignment
from read_reference import Read_reference
from genome_reference import Genome_reference
from complement import Complement

def readargs(args=None):
    parser = argparse.ArgumentParser(prog='Biosequence project',)
    #General arguments
    parser.add_argument('-t', '--task',help="task",required=True)
    parser.add_argument('-g', '--genomefile',help="Genome fasta file (multiple records)",)
    parser.add_argument('-r', '--referencefile',help="kdb file. Can be either input or name for output file",)
    parser.add_argument('-k', '--kmer-size',help="length of kmers",)
    #Task specific arguments
    #align
    parser.add_argument('-a', '--alignfile',help="aln file. Can be either input or name for output file",)
    parser.add_argument('--reads',help="fastq reads file",)
    parser.add_argument('-m', '--unique-threshold',help="unique k-mer threshold",)
    parser.add_argument('-p', '--ambiguous-threhold',help="ambiguous k-mer threshold",)               
    #align+rc
    parser.add_argument('--reverse-complement',action='store_true',)
    #align+quality
    parser.add_argument('--min-read-quality',type=int)
    parser.add_argument('--min-kmer-quality',type=int,)
    parser.add_argument('--max-genomes',type=int,)
    #coverage
    parser.add_argument('--genomes',)
    parser.add_argument('--regions',)
    parser.add_argument('--window-size',type=int,default=100,)
    parser.add_argument('--min-coverage',type=int,default=1,)
    parser.add_argument('--full-coverage',action='store_true',)
    #similarity
    parser.add_argument('--filter-similar',action='store_true',)
    parser.add_argument('--similarity-threshold',type=float,default=0.95,)
    
    return parser.parse_args(args)
 


def save_reference_file(data:Any, filename:str)-> None:
    with open(filename, 'wb') as file:
        pickle.dump(data, file, protocol=pickle.HIGHEST_PROTOCOL)
    print("file saved as", filename)

def load_reference_file(filename: str):
    with open(filename, 'rb') as file:
        return pickle.load(file)

def build_genome_reference(genome_file:str, k:int , ref_file: Optional[str], filter_similar: Optional[bool], similarity_threshold: Optional[float])-> Optional[Genome_reference]:
    genome_ref = Genome_reference(genome_file,k)
    if filter_similar:
        if not similarity_threshold:
            similarity_threshold = 0.95
        genome_ref.filter_by_similarity(similarity_threshold)
    if ref_file :
        save_reference_file(genome_ref, ref_file)
        print(f"Reference file saved as {ref_file}")
    else:
         return genome_ref

def dump_genome(genome_ref:Genome_reference, filter_similar: Optional[bool])->None:
    total_kmer_dict ={}
    genome_dict = {}
    for kmer, genome in genome_ref.create_kmer_ref().items():
            total_kmer_dict[kmer]= genome
    genome_ref.kmer_per_genome()
    for genome in genome_ref.genome_dictionary.items():
        genome_dict[genome[0]]= {"total_bases": genome[1].length, "unique_kmers": genome[1].specific, "multi_mapping_kmers": genome[1].multimapping}
    if filter_similar:
        dump_dict = {"Kmers": total_kmer_dict, "Summary": genome_dict},genome_ref.similarity_dict
    else:
        dump_dict = {"Kmers": total_kmer_dict, "Summary": genome_dict}
    print(json.dumps(dump_dict, indent=1))

def align(reads:str, genome_ref:Genome_reference, p:int ,m:int, align_file: Optional[str], min_read_quality: Optional[int] ,min_kmer_quality: Optional[int], max_genome_num: Optional[int], comp:str ="regular"):
    skipped = [0,0,0]
    read_ref = Read_reference(reads,genome_ref.k)
    if min_read_quality:
        skipped[0] = read_ref.filter_by_read_quality(min_read_quality)
    if comp == "complement":
        match = Matching_algorythme( read_ref.reads_dict_to_array_complement(), read_ref.read_dictionary, genome_ref.genome_dict_to_lst_of_array(), genome_ref.genome_dictionary)
    else:
        match = Matching_algorythme( read_ref.read_array, read_ref.read_dictionary, genome_ref.genome_dict_to_lst_of_array(), genome_ref.genome_dictionary)
    align = Pseudo_alignment(match.panda_data_frame, read_ref.read_number, read_ref.read_dictionary ,int(p), int(m))
    if  min_kmer_quality:
        skipped[1] = align.filter_by_read_quality_of_kmer(min_kmer_quality,genome_ref.k)
    if max_genome_num:
        skipped[2] = align.filter_by_read_apearance(max_genome_num)
    if max_genome_num or min_read_quality or min_kmer_quality:
        align.skip_setter(skipped)
    if align_file:
        save_reference_file(align, align_file)
        return None
    else:
        return align

def dump_align(align,min_read_quality ,min_kmer_quality ,max_genome_num,reverse_complement)->None:
    #if reverse complement using the regular algo for the winner while implementing the looser in the summary per genome
    if reverse_complement:
        comp:Complement = align
        align = comp.winner
    spe = align.spe_len
    multimapping = align.multi_len
    total_read = align.read_number
    genome_dict = {}
    #making a double genomedict for the revers complement
    if reverse_complement:
        for genome in align.aligne():
            genome_dict[genome+"_F"]= {"unique_reads":align.specifique_df[align.specifique_df["Genome"] == genome]["Read"].nunique(),"ambiguous_reads":align.multimapping_df[align.multimapping_df["Genome"] == genome]["Read"].nunique()}
            genome_dict[genome+"_R"]= {"unique_reads":comp.looser.specifique_df[comp.looser.specifique_df["Genome"] == genome]["Read"].nunique(),"ambiguous_reads":comp.looser.multimapping_df[comp.looser.multimapping_df["Genome"] == genome]["Read"].nunique()}
    else:
        for genome in align.aligne():
            genome_dict[genome]= {"unique_reads":align.specifique_df[align.specifique_df["Genome"] == genome]["Read"].nunique(),"ambiguous_reads":align.multimapping_df[align.multimapping_df["Genome"] == genome]["Read"].nunique()}
    #addapting the Statistics in case the one of the quality filter was applied
    if min_read_quality or min_kmer_quality or max_genome_num:
        dump_dict = {"Statistics":{"unique_mapped_reads":spe ,"ambiguous_mapped_reads":multimapping,"unmapped_reads":total_read-(spe+multimapping-align.skip[0]-align.skip[1]-align.skip[2]),"filtered_quality_reads":align.skip[0],"filtered_quality_kmers": align.skip[1],"filtered_hr_kmers": align.skip[2]},"Summary": genome_dict}
    else:
        dump_dict = {"Statistics":{"unique_mapped_reads":spe ,"ambiguous_mapped_reads":multimapping,"unmapped_reads":total_read-(spe+multimapping)},"Summary": genome_dict}
    print(json.dumps(dump_dict, indent=1))

    
def main():
    args = readargs()
    
    if args.task == "reference":
        if args.genomefile and args.kmer_size and args.referencefile:
            build_genome_reference(args.genomefile, args.kmer_size, args.referencefile,args.filter_similar,args.similarity_threshold)  
        
    if args.task == "dumpref":
        if args.genomefile and args.kmer_size:
            gen_ref = build_genome_reference(args.genomefile, args.kmer_size, args.referencefile,args.filter_similar,args.similarity_threshold)
            dump_genome(gen_ref,args.filter_similar)
        elif args.referencefile:
             dump_genome(load_reference_file(args.referencefile),args.filter_similar)

    if args.task == "align":
        if args.reads and args.alignfile:
            if not args.unique_threshold:
                args.unique_threshold = 1
            if not args.ambiguous_threhold:
                args.ambiguous_threhold = 1
            if args.genomefile and args.kmer_size:
                genome_ref = Genome_reference(args.genomefile, args.kmer_size)
            elif args.referencefile:
                genome_ref :Genome_reference= load_reference_file(args.referencefile)
            if args.reverse_complement:
                regular = align(args.reads, genome_ref, args.ambiguous_threhold , args.unique_threshold, None,args.min_read_quality,args.min_kmer_quality,args.max_genomes,"regular")
                complement = align(args.reads, genome_ref, args.ambiguous_threhold , args.unique_threshold, None,args.min_read_quality,args.min_kmer_quality,args.max_genomes,"complement")
                comp = Complement(regular,complement)
                save_reference_file(comp,args.alignfile)
            else:
                align(args.reads, genome_ref, args.ambiguous_threhold , args.unique_threshold, args.alignfile,args.min_read_quality,args.min_kmer_quality,args.max_genomes)
    if args.task == "dumpalign":
        if args.genomefile and args.kmer_size and args.reads:
            genome_ref = build_genome_reference(args.genomefile, args.kmer_size,args.referencefile,args.filter_similar,args.similarity_threshold)
        if args.reads and args.referencefile :
            genome_ref= load_reference_file(args.referencefile)
        if not args.alignfile:
            if not args.unique_threshold:
                args.unique_threshold = 1
            if not args.ambiguous_threhold:
                args.ambiguous_threhold = 1
            if args.reverse_complement:
                regular = align(args.reads, genome_ref, args.ambiguous_threhold , args.unique_threshold, None,args.min_read_quality,args.min_kmer_quality,args.max_genomes,"regular")
                complement = align(args.reads, genome_ref, args.ambiguous_threhold , args.unique_threshold, None,args.min_read_quality,args.min_kmer_quality,args.max_genomes,"complement")
                comp = Complement(regular,complement)
            else:
                aligne = align(args.reads, genome_ref, args.ambiguous_threhold , args.unique_threshold, args.alignfile,args.min_read_quality,args.min_kmer_quality,args.max_genomes)
        else:
            if args.reverse_complement:
                regular = align(args.reads, genome_ref, args.ambiguous_threhold , args.unique_threshold, None,args.min_read_quality,args.min_kmer_quality,args.max_genomes,"regular")
                complement = align(args.reads, genome_ref, args.ambiguous_threhold , args.unique_threshold, None,args.min_read_quality,args.min_kmer_quality,args.max_genomes,"complement")
                comp = Complement(regular,complement)
            else:
                aligne = load_reference_file(args.alignfile)
        if args.reverse_complement:
            dump_align(comp,args.min_read_quality,args.min_kmer_quality,args.max_genomes,args.reverse_complement)
        else:
            dump_align(aligne,args.min_read_quality,args.min_kmer_quality,args.max_genomes,args.reverse_complement)
    
    
    

if __name__=="__main__":
    main()

