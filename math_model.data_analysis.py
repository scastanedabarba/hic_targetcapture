import pysam
import pandas as pd
import glob
import os

def extract_links(bam_file):
    """
    Function to be applied to file after matlock filtering, no unaligned reads should be present.
    Takes a bam file as input and returns the counts for the links present.
    """
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    count = 0
    genomes = {'pB10::rfp_putative':'pB10', 'NC_002947.4':'P.Putida', 'NC_000913.3_gfp':'E.Coli'}
    pB10_pB10 = 0
    pB10_ecoli = 0
    pB10_pputida = 0
    ecoli_ecoli = 0
    pputida_pputida = 0
    multimapped = 0
    pairs = []
    
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
                genome='clipped'
        pairs.append([name, genome, loc, length, read.get_tags(), read])
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
                if links == 'pB10P.Putida':
                    pB10_pputida += 1
                if links == 'P.PutidapB10':
                    pB10_pputida += 1
                if links == 'pB10E.Coli':
                    pB10_ecoli += 1
                if links == 'E.ColipB10':
                    pB10_ecoli += 1
                if links == 'E.ColiE.Coli':
                    ecoli_ecoli += 1 
                if links == 'P.PutidaP.Putida':
                    pputida_pputida += 1
                pairs = []
    return([pB10_pputida, pB10_ecoli, ecoli_ecoli, pB10_pB10, pputida_pputida])

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


os.chdir('/mnt/ceph/cast9836/00_projects/hic_targetcapture/05_bwa_alignments/math_model/sorted/')
files = glob.glob('*.bam')
files.sort()

lib = 0
total = len(files)

df = pd.DataFrame(columns = ['Library', 'Depth', 'Trial', 'pB10_pputida'])

for file in files:
    lib += 1
    index = file[0]
    depth = file[2:10]
    trial = file[11]
    
    counts = extract_links(file)
    df = df.append({'Library': index, 'Depth':depth, 'Trial':trial, 'pB10_pputida':counts[0]}, ignore_index=True)
    
    percent = (lib/total)*100            
    print('Library ', index, ' at ', depth, ' sequencing depth and trial ', trial, ' done ::: ', round(percent,2))
    
    
df.to_csv("/mnt/ceph/cast9836/00_projects/hic_targetcapture/06_output/math_model/read_simulations.csv", index=False)