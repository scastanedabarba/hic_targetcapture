import pysam
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import seaborn as sns
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker
import copy
import random
from Bio import SeqIO
from hic_functions import extract_links, update_counts, check_alt

def extract_links(bam_file):
    """

    """
    print('Starting ', bam_file)
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    count = 0
    genomes = {'pB10::rfp_putative':'pB10', 'NC_002947.4':'P.Putida', 'NC_000913.3_gfp':'E.Coli'}
    pB10_pB10 = 0
    pB10_ecoli = 0
    pB10_pputida = 0
    ecoli_ecoli = 0
    pputida_pputida = 0
    pairs = []
    
    #lists that store the location on the pB10 genome of each aligned read
    pB10_pB10_loc = []
    pB10_pputida_loc = []
    pB10_ecoli_loc = []
    
    #lists that store location on chromosome of pb10 aligned reads
    ecoli_pB10_loc = []
    pputida_pB10_loc = []
    
    global p_putida
    global p_p
    global p_ecoli
    p_putida = [0]*pb10_len
    p_p = [0]*pb10_len
    p_ecoli = [0]*pb10_len
    
    for read in samfile.fetch(until_eof=True):
        count += 1
        genome = ''
        name = ''
        cigar = ''
        loc = ''
        seq = read.query_sequence
        length = len(seq)
        if read.is_unmapped:
            name = read.query_name
            genome = 'unaligned'
        if not read.is_unmapped:
            name = read.query_name
            genome = genomes.get(read.reference_name)
            cigar = read.cigarstring
            match = str(length)+'M'
            loc = read.reference_start
            if cigar != match:
                genome='unaligned'
        pairs.append([name, genome, loc, length, read.get_tags(), read.query_sequence])
        if count % 2 == 0:
            r1_multimapped = check_alt(pairs[0][1], pairs[0][4], pairs[0][3])
            r2_multimapped = check_alt(pairs[1][1], pairs[1][4], pairs[1][3])
            if len(r1_multimapped) > 1 or len(r2_multimapped) >1:
                multimapped += 1
                pairs = []
            else:
                links = pairs[0][1] + pairs[1][1]
                if links == 'pB10pB10':
                    pB10_pB10 += 1
                    update_counts('plasmid','plasmid', [pairs[0][2], pairs[1][2]], [pairs[0][3], pairs[1][3]])
                    pB10_pB10_loc.append(pairs[0][2])
                    pB10_pB10_loc.append(pairs[1][2])
                if links == 'pB10P.Putida':
                    pB10_pputida += 1
                    update_counts('plasmid', 'putida', pairs[0][2], pairs[0][3])
                    pB10_pputida_loc.append(pairs[0][2])
                    pputida_pB10_loc.append(pairs[1][2])
                    #pputida_file.write(">" + pairs[1][0] + '\n' + pairs[1][-1] + '\n')
                if links == 'P.PutidapB10':
                    pB10_pputida += 1
                    update_counts('putida', 'plasmid', pairs[1][2], pairs[1][3])
                    pB10_pputida_loc.append(pairs[1][2])
                    pputida_pB10_loc.append(pairs[0][2])
                    #pputida_file.write(">" + pairs[0][0] + '\n' + pairs[0][-1] + '\n')
                if links == 'pB10E.Coli':
                    pB10_ecoli += 1
                    update_counts('plasmid', 'ecoli', pairs[0][2], pairs[0][3])
                    pB10_ecoli_loc.append(pairs[0][2])
                    ecoli_pB10_loc.append(pairs[1][2])
                    #ecoli_file.write(">" + pairs[1][0] + '\n' + pairs[1][-1] + '\n')
                if links == 'E.ColipB10':
                    pB10_ecoli += 1
                    update_counts('ecoli', 'plasmid', pairs[1][2], pairs[1][3])
                    pB10_ecoli_loc.append(pairs[1][2])
                    ecoli_pB10_loc.append(pairs[1][2])
                    #ecoli_file.write(">" + pairs[0][0] + '\n' + pairs[0][-1] + '\n')
                if links == 'E.ColiE.Coli':
                    ecoli_ecoli += 1 
                if links == 'P.PutidaP.Putida':
                    pputida_pputida += 1
                pairs = []
    print(np.array([count/2, multimapped, pB10_pB10, pB10_ecoli, pB10_pputida], '\n'))
    print('\n')
    return(np.array([[count/2, multimapped, pB10_pB10, pB10_ecoli, pB10_pputida],
                    [p_putida, p_p, p_ecoli], 
                    [pB10_pputida_loc, pB10_pB10_loc, pB10_ecoli_loc, ecoli_pB10_loc, pputida_pB10_loc]]))

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
        if link == 'putida':
            s = start_index
            e = start_index + length
            for i in range(s, e):
                if i < pb10_len:
                    p_putida[i] += 1
                else:
                    i = i - pb10_len
                    p_putida[i] += 1 
        

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
        if ref == 'putida':
            s = start_index
            e = start_index + length
            for i in range(s, e):
                if i < pb10_len:
                    p_putida[i] += 1
                else:
                    i = i - pb10_len
                    p_putida[i] += 1 
                    
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



os.chdir('/mnt/ceph/cast9836/00_projects/hic_targetcapture/05_bwa_alignment/filtered/')

pb10_len = 68345
p_p = [0]*pb10_len
p_ecoli = [0]*pb10_len
p_putida = [0]*pb10_len

A_all = extract_links('B.bam')
B_all = extract_links('C.bam')
C_all = extract_links('D.bam')
D_all = extract_links('E.bam')
E_all = extract_links('F.bam')
ctrl1_all = extract_links('ctrl1_ecoli.bam')

Bplus_all = extract_links('C+.bam')
Cplus_all = extract_links('D+.bam')
Dplus_all = extract_links('E+.bam')
Eplus_all = extract_links('F+.bam')
Fplus_all = extract_links('G+.bam')
ctrl1plus_all = extract_links('ctrl1+_ecoli.bam')
ctrl2plus_all = extract_links('ctrl2+_soil.bam')

#get sequencing depth for each library from multiqc report
os.chdir('/mnt/ceph/cast9836/00_projects/hic_targetcapture/02_fastqc/untrimmed/multiqc_data/')
seq_depth = pd.read_csv('multiqc_general_stats.txt', sep='\t')
seq_depth = seq_depth[['Sample', 'FastQC_mqc-generalstats-fastqc-total_sequences']]
seq_depth.rename(columns={'FastQC_mqc-generalstats-fastqc-total_sequences':'num_reads'}, inplace=True)
seq_depth = seq_depth.iloc[::2]
seq_depth['Sample'] = seq_depth.Sample.str[:-3]
seq_depth.reset_index(inplace=True)
seq_depth.drop(columns=['index'], inplace=True)
seq_depth['Type'] = ['Hi-C', 'Hi-C+', 'Hi-C', 'Hi-C+', 'Hi-C', 'Hi-C+', 'Hi-C', 'Hi-C+', 'Hi-C', 'Hi-C+', 'Hi-C+', 'Hi-C', 'Hi-C+']
seq_depth['Sample'] = ['A', 'B', 'B', 'C', 'C', 'D', 'D', 'E', 'E', 'F', 'Ctrl1', 'Ctrl1', 'Ctrl2']
seq_depth


#place hic-data into dataframe, output two csv files
#first is raw counts and second is normalized per million counts
hic_all = pd.DataFrame([A_all[0], B_all[0], C_all[0], D_all[0], E_all[0], ctrl1_all[0]],
                     columns = ['count', 'multimapped', 'pB10_pB10','pB10_ecoli','pB10_pputida'])
hic_all['Sample'] = ['A', 'B', 'C', 'D', 'E', 'Ctrl1']
hic_all = hic_all.merge(seq_depth[seq_depth['Type']=='Hi-C']).drop(['count', 'Type'], axis=1)
hic_all.set_index("Sample", inplace=True)
hic_all.to_csv('/mnt/ceph/cast9836/00_projects/hic_targetcapture/06_output/hic_counts_raw.csv')

hic_all['standard'] = seq_depth['num_reads'].min()/hic_all['num_reads']
hic_all.drop(columns='num_reads', inplace=True)

#normalize hic counts
hic_all['pB10_pB10'] = hic_all['pB10_pB10']*hic_all['standard']
hic_all['pB10_ecoli'] = hic_all['pB10_ecoli']*hic_all['standard']
hic_all['pB10_pputida'] = hic_all['pB10_pputida']*hic_all['standard']

hic_all.drop(columns='standard', inplace=True)
hic_all.to_csv('/mnt/ceph/cast9836/00_projects/hic_targetcapture/06_output/hic_counts_normalized.csv')

#place hic+-data into dataframe, output two csv files
#first is raw counts and second is normalized per million counts
hicplus_all = pd.DataFrame([Bplus_all[0], Cplus_all[0], Dplus_all[0], Eplus_all[0], Fplus_all[0], ctrl1plus_all[0], ctrl2plus_all[0]],
                     columns = ['count', 'multimapped', 'pB10_pB10','pB10_ecoli','pB10_pputida'])
hicplus_all['Sample'] = ['B', 'C', 'D', 'E', 'F', 'Ctrl1', 'Ctrl2']
hicplus_all = hicplus_all.merge(seq_depth[seq_depth['Type']=='Hi-C+']).drop(['count', 'Type'], axis=1)
hicplus_all.set_index("Sample", inplace=True)
hicplus_all.to_csv('/mnt/ceph/cast9836/00_projects/hic_targetcapture/06_output/hic+_counts_raw.csv')

hicplus_all['standard'] = seq_depth['num_reads'].min()/hicplus_all['num_reads']
norm_factor = list(hicplus_all['standard'])
print(hicplus_all)
hicplus_all.drop(columns='num_reads', inplace=True)

#normalize hic counts
hicplus_all['pB10_pB10'] = hicplus_all['pB10_pB10']*hicplus_all['standard']
hicplus_all['pB10_ecoli'] = hicplus_all['pB10_ecoli']*hicplus_all['standard']
hicplus_all['pB10_pputida'] = hicplus_all['pB10_pputida']*hicplus_all['standard']

hicplus_all.drop(columns='standard', inplace=True)
hicplus_all.to_csv('/mnt/ceph/cast9836/00_projects/hic_targetcapture/06_output/hic+_counts_normalized.csv')