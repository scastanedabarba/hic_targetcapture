import pysam
import pandas as pd
import numpy as np
import sys

def process_reads(bam_file, out_name):
    
    samfile = pysam.AlignmentFile(bam_file, 'rb')

    genomes = {'pB10::rfp_putative':'pB10', 'NC_000913.3_gfp':'E.Coli'}
    pairs = []
    count=0
        
    data = []
    
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
        pairs.append([name, genome, loc, length, read.get_tags(), seq])
        if count % 2 == 0:
            #first, check for multimapped reads and remove pairs where one read is multimapped
            r1_multimapped = check_alt(pairs[0][1], pairs[0][4], pairs[0][3])
            r2_multimapped = check_alt(pairs[1][1], pairs[1][4], pairs[1][3])
            if len(r1_multimapped) > 1 or len(r2_multimapped) >1:
                pairs = []
            else:
                link = pairs[0][1]+pairs[1][1]
                if link == 'unalignedunaligned': #don't append unaligned-unaligned reads to df
                    pass
                else:
                    data.append([pairs[0][0], link, pairs[0][5], pairs[0][3], 
                                 pairs[0][2], pairs[1][5], pairs[1][3], pairs[1][2]])
            pairs=[]
    data_df = pd.DataFrame(data, columns=['seqid', 'link', 'R1','R1length','R1loc','R2','R2length','R2loc'])
    data_df.to_csv(out_name)
    return(data_df)

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
        
#assign command line input to variables
file = sys.argv[1]
out_df = sys.argv[2]
out_pothosts = sys.argv[3]

#create df with all reads and their alignment
all_reads = process_reads(file, out_df)

#output potential pB10 hosts to 
extract_pothosts(all_reads, out_pothosts)
                