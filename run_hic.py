def extract_links(bam_file):
    """
    Function to be applied to file after matlock filtering, no unaligned reads should be present.
    Takes a bam file as input and returns the counts for the links present.
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
        if count % 1000000 == 0:
            print(np.array([count/2, multimapped, pB10_pB10, pB10_ecoli, pB10_pputida]))
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
        pairs.append([name, genome, loc, length, read.get_tags(), read])
        if count % 2 == 0:
            r1_multimapped = check_alt(pairs[0][1], pairs[0][4], pairs[0][3])
            r2_multimapped = check_alt(pairs[1][1], pairs[1][4], pairs[1][3])
            if len(r1_multimapped) > 1 or len(r2_multimapped) >1:
                #print(pairs[0][1] + pairs[1][1])
                #print(pairs[0][-1],'\n',pairs[1][-1])
                multimapped += 1
            else:
                links = pairs[0][1] + pairs[1][1]
                if links == 'pB10pB10':
                    pB10_pB10 += 1
                    #update_counts('plasmid','plasmid', [pairs[0][2], pairs[1][2]], [pairs[0][3], pairs[1][3]])
                if links == 'pB10P.Putida':
                    pB10_pputida += 1
                    #update_counts('plasmid', 'putida', pairs[0][2], pairs[0][3])
                if links == 'P.PutidapB10':
                    pB10_pputida += 1
                    #update_counts('putida', 'plasmid', pairs[1][2], pairs[1][3])
                if links == 'pB10E.Coli':
                    pB10_ecoli += 1
                    #update_counts('plasmid', 'ecoli', pairs[0][2], pairs[0][3])
                if links == 'E.ColipB10':
                    pB10_ecoli += 1
                    #update_counts('ecoli', 'plasmid', pairs[1][2], pairs[1][3])
                if links == 'E.ColiE.Coli':
                    ecoli_ecoli += 1 
                if links == 'P.PutidaP.Putida':
                    pputida_pputida += 1
                pairs = []
    print(np.array([count/2, multimapped, pB10_pB10, pB10_ecoli, pB10_pputida], '\n'))
    print('\n')
    return(np.array([count/2, pB10_pB10, pB10_ecoli, pB10_pputida]))


def update_counts(ref, link, start_index, length):
    """Function for updating the coverage list"""
    global p_unaligned
    global p_ecoli
    global p_p
    
    if ref == 'plasmid':
        if link == 'unaligned':
            s = start_index
            e = start_index + length
            for i in range(s, e):
                if i < pb10_len:
                    p_unaligned[i] += 1
                else:
                    i = i - pb10_len
                    p_unaligned[i] += 1
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
        if ref == 'unaligned':
            s = start_index
            e = start_index + length
            for i in range(s, e):
                if i < pb10_len:
                    p_unaligned[i] += 1
                else:
                    i = i - pb10_len
                    p_unaligned[i] += 1
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

