## Originally written by Justin Tan and modified by Patrick Phaneuf

import pysam as ps
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from matplotlib import pyplot as plt
from scipy import stats, optimize


def geneFrame(genbank_file):
    """Function: Parses genbank file to dataframe format for easier downstream utilization
	Input: Genbank file
    Output: Dataframe containing all genes, start/stop locations and strand (plus/minus), 
    as well as gene function, product and amino acid sequence
    """
    from Bio import SeqIO
    
    infile = SeqIO.read(genbank_file,'gb')
    genes = []
    name = []
    product = []
    func = []
    strand = []
    start = []
    stop = []
    aaseq = []
    cds_seq = []
    codon_start = []
    
    genome_seq_df = pd.DataFrame({'sequence':list(str(infile.seq))},index=range(1,len(str(infile.seq))+1))
    
    for feature in infile.features:
        
        if feature.type == 'CDS' and 'product' in feature.qualifiers:  
            #Only cares for coding sequences which are not pseudogenes
            
            try: genes.append(feature.qualifiers['locus_tag'][0])
            except: genes.append('')
            
            try: name.append(feature.qualifiers['gene'][0])
            except: name.append('')

            try: codon_start.append(int(feature.qualifiers['codon_start'][0]))
            except: codon_start.append(np.nan)
                
            product.append(feature.qualifiers['product'][0])
            cds_seq.append(str(feature.location.extract(infile.seq)))
            
            if 'function' in feature.qualifiers: #not all genes have known functions
                func.append(feature.qualifiers['function'][0])
            else:
                func.append("N/A")
                
            try: aaseq.append(feature.qualifiers['translation'][0])
            except: aaseq.append("N/A")
                
            if feature.strand == 1:
                strand.append("plus") 
                start.append(feature.location.start.real+1)  
                stop.append(feature.location.end.real)
            elif feature.strand == -1:
                strand.append("minus") 
                start.append(feature.location.start.real+1)
                stop.append(feature.location.end.real)
                
    df = pd.DataFrame({"gene": genes, "name": name, "product": product, "function": func, "strand": strand, "start": start, "stop": stop, "codon_start":codon_start, "cds_seq":cds_seq,"aaseq": aaseq},
                          columns = ["gene", "name", "function", "product", "strand", "start", "stop", "codon_start","cds_seq","aaseq"])
    df = df.set_index("gene")
    df['id'] = (infile.id).split('.')[0]
    return df, genome_seq_df
