import pysam
import pandas as pd
import numpy as np
import sys
import os

def extract_links(bam_file, out_name):
    """
    Function to be applied to file after matlock filtering, no unaligned reads should be present.
    Takes a bam file as input and returns the counts for the links present.
    """
    print('Starting ', bam_file)
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    count = 0
    genomes = {'pB10::rfp_putative':'pB10', 'NC_000913.3_gfp':'E.Coli'}
    pB10_pB10 = 0
    pB10_ecoli = 0
    pB10_unaligned = 0
    pairs = []
    
    #lists that store the location on the pB10 genome of each aligned read
    pB10_pB10_loc = []
    pB10_ecoli_loc = []
    
    #lists that store location on chromosome of pb10 aligned reads
    ecoli_pB10_loc = []
    
    data = []
    
    global p_unaligned
    global p_p
    global p_ecoli
    p_unaligned = [0]*pb10_len
    p_p = [0]*pb10_len
    p_ecoli = [0]*pb10_len
    
    for read in samfile.fetch(until_eof=True):
        count += 1
        genome = ''
        name = ''
        #cigar = ''
        loc = ''
        seq = read.query_sequence
        length = len(seq)
        if read.is_unmapped:
            name = read.query_name
            genome = 'unaligned'
        if not read.is_unmapped:
            name = read.query_name
            genome = genomes.get(read.reference_name)
            #cigar = read.cigarstring
            loc = read.reference_start
            #match = str(length)+'M'
            #if cigar != match:
            #    genome='unaligned'
        pairs.append([name, genome, loc, length, read.get_tags(), read.query_sequence])
        if count % 2 == 0:
            r1_multimapped = check_alt(pairs[0][1], pairs[0][4], pairs[0][3])
            r2_multimapped = check_alt(pairs[1][1], pairs[1][4], pairs[1][3])
            if len(r1_multimapped) > 1 or len(r2_multimapped) >1:
                pairs = []
            else:
                links = pairs[0][1] + pairs[1][1]
                if links == 'pB10pB10':
                    pB10_pB10 += 1
                    update_counts('plasmid','plasmid', [pairs[0][2], pairs[1][2]], [pairs[0][3], pairs[1][3]])
                    pB10_pB10_loc.append(pairs[0][2])
                    pB10_pB10_loc.append(pairs[1][2])
                if links == 'pB10E.Coli':
                    pB10_ecoli += 1
                    update_counts('plasmid', 'ecoli', pairs[0][2], pairs[0][3])
                    pB10_ecoli_loc.append(pairs[0][2])
                    ecoli_pB10_loc.append(pairs[1][2])
                if links == 'E.ColipB10':
                    pB10_ecoli += 1
                    update_counts('ecoli', 'plasmid', pairs[1][2], pairs[1][3])
                    pB10_ecoli_loc.append(pairs[1][2])
                    ecoli_pB10_loc.append(pairs[1][2])
                if links == 'pB10unaligned':
                    pB10_unaligned += 1
                    data.append([pairs[0][0], links, pairs[0][5], pairs[0][3], 
                                 pairs[0][2], pairs[1][5], pairs[1][3], pairs[1][2]])
                if links == 'unalignedpB10':
                    pB10_unaligned += 1
                    data.append([pairs[0][0], links, pairs[0][5], pairs[0][3], 
                                 pairs[0][2], pairs[1][5], pairs[1][3], pairs[1][2]])
                pairs = []                   
                
    data_df = pd.DataFrame(data, columns=['seqid', 'link', 'R1','R1length','R1loc','R2','R2length','R2loc'])
    data_df.to_csv(out_name+'_pb10_unaligned.csv')
    return([[count/2, pB10_pB10, pB10_ecoli, pB10_unaligned], 
            [p_p, p_ecoli], 
            [pB10_pB10_loc, pB10_ecoli_loc, ecoli_pB10_loc]], data_df)

def update_counts(ref, link, start_index, length):
    """Function for updating the coverage list"""
    global p_putida
    global p_ecoli
    global p_p
    
    if ref == 'plasmid':
        if link == 'ecoli':
            s = start_index
            e = start_index + length
            for i in range(s, e):
                if i < pb10_len:
                    p_ecoli[i] += 1
                else:
                    i = i - pb10_len
                    p_ecoli[i] += 1         

    if link == 'plasmid':
        if ref == 'ecoli':
            s = start_index
            e = start_index + length
            for i in range(s, e):
                if i < pb10_len:
                    p_ecoli[i] += 1
                else:
                    i = i - pb10_len
                    p_ecoli[i] += 1 
                    
    if ref == 'plasmid' and link == 'plasmid':                  
        for i in range(0,2):
            s = start_index[i]
            e = start_index[i] + length[i]
            for i in range(s, e):
                if i < pb10_len:
                    p_p[i] += 1
                else:
                    i = i - pb10_len
                    p_p[i] += 1
                    
def check_alt(aligned, tags, length):
    """
    This function uses the alignment tags to determine if a read was
    mapped perfectly to multiple references. The output is a list of 
    references to which the read mapped.
    """
    genomes = [aligned]
    for tag in tags:
        if tag[0] == 'XA' or tag[0] == 'SA':
            tag = list(tag)
            others = tag[1].split(';')[:-1]
            fmatch = []
            for alt in others:
                info = alt.split(',')
                if info[-2] == str(length)+'M':
                    fmatch.append(info[0])
            for i in fmatch:
                genomes.append(i)
    return(list(set(genomes)))

def extract_pothosts(df, out_name):
    output_file = open(out_name,'w')
    set1_seqid = df[df['link']=='pB10unaligned']['seqid'].tolist()
    set1_sequence = df[df['link']=='pB10unaligned']['R2'].tolist()
    set2_seqid = df[df['link']=='unalignedpB10']['seqid'].tolist()
    set2_sequence = df[df['link']=='unalignedpB10']['R1'].tolist()
    all_seqid = set1_seqid + set2_seqid
    all_sequence = set1_sequence + set2_sequence
    for seqid, seq in zip (all_seqid, all_sequence):
        output_file.write(">" + seqid + '\n' + seq + '\n') 

pb10_len = 68345
p_p = [0]*pb10_len
p_ecoli = [0]*pb10_len
p_unaligned = [0]*pb10_len

os.chdir('/mnt/ceph/cast9836/00_projects/hic_targetcapture/05_bwa_alignments/denovo/sorted/')
input_file = sys.argv[1]
file = input_file.split('_')[0]

results, reads = extract_links(str(input_file), file)

os.chdir('/mnt/ceph/cast9836/00_projects/hic_targetcapture/06_output/denovo/')

read_count = pd.DataFrame({'pb10_pb10': [results[0][1]],
                           'pb10_ecoli': [results[0][2]],
                           'pB10_unaligned': [results[0][3]]})
read_count.to_csv(file+'_read_count.csv', index=False)

plasmid_coverage = pd.DataFrame({'pb10_pb10':results[1][0], 
                                 'pb10_ecoli':results[1][1]})
plasmid_coverage.to_csv(file+'_plasmid_coverage.csv', index=False)


with open(file+'_locations.txt', 'w') as file:
    file.write('pB10_pB10_loc, pB10_ecoli_loc, ecoli_pB10_loc\n')
    for item in results[2]:
        file.write(str(item))
        file.write('\n')

extract_pothosts(reads, file+"_pothosts.fasta")