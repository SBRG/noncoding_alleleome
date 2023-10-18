NON_K12_EXP_L = ['LTEE', 'LTEE_Ara']


# Needs to reflect the opposite transformations into the BOP27 genome.
def get_BOP27_pos_from_K12_pos(K12_pos):
    BOP27_pos = K12_pos
    # apply each of the following transformations according to position
    
    # Deletion of 2 basepairs at 2173363 (gatC)
    if K12_pos == 2173363 or K12_pos == 2173364:
        print(str(K12_pos) + " position no longer exists in the BOP27 genome.")
        BOP27_pos = -1
    if K12_pos >= 2173365:
        BOP27_pos -= 2
    if K12_pos >= 3560456:
        BOP27_pos += 1
    if K12_pos >= 4296381:
        BOP27_pos += 2
        
    return BOP27_pos


def get_K12_pos_from_BOP27(BOP27_pos):
    K12_pos = BOP27_pos
    
    if BOP27_pos >= 2173365:
        K12_pos += 2
    if BOP27_pos >= 3560456:
        K12_pos -= 1
    if BOP27_pos >= 4296381:
        K12_pos -= 2
        
    return K12_pos


def is_overlap(range1, range2):
    is_overlap = False
    if range1 and range2:
        r1 = range(range1[0], range1[1] + 1)  # Python's range function doesn't include final number given to it.
        r2 = range(range2[0], range2[1] + 1)  # Python's range function doesn't include final number given to it.
        is_overlap = bool(set(r1).intersection(r2))
    return is_overlap


# used to be get_site_hit_set
def get_feature_hit_set(mut_range, struct_df, range_col, object_hit_col):
    struct_df_copy = struct_df.copy()
    struct_df_copy["mutation hit"] = struct_df_copy[range_col].apply(is_overlap, args=(mut_range, ))
    hit_df = struct_df_copy[struct_df_copy["mutation hit"]]
    hit_set = set(hit_df[object_hit_col])
    return hit_set


# http://regulondb.ccg.unam.mx/menu/using_regulondb/tutorials/project_glossary/index.jsp
# A promoter is the DNA sequence where RNA polymerase binds and initiates transcription.
# Notes: Promoter sequences are specific to the different sigma factors associated to the RNA polymerase core.
# A promoter is represented as a stretch of 60 upstream and 20 downstream upper-case nucleotide sequences from the precise initiation of transcription, also called +1.
REGULONDB_PROMOTER_TSS_DOWNSTREAM_LEN = 20
REGULONDB_PROMOTER_TSS_UPSTREAM_LEN = 60


import pandas as pd


# Assuming genome positions are 1 based.
def get_promoter_range_from_RegulonDB_df_row(promoter_df_row):
	r = ()
	promoter_seq_str = promoter_df_row[9]
	promoter_strand = promoter_df_row[2]
	if not pd.isna(promoter_seq_str):
		TSS_genome_pos = int(promoter_df_row[3])
		r = (TSS_genome_pos - REGULONDB_PROMOTER_TSS_UPSTREAM_LEN, TSS_genome_pos + REGULONDB_PROMOTER_TSS_DOWNSTREAM_LEN)
		if promoter_strand == "reverse":
			r = (TSS_genome_pos - REGULONDB_PROMOTER_TSS_DOWNSTREAM_LEN, TSS_genome_pos + REGULONDB_PROMOTER_TSS_UPSTREAM_LEN)
	return r
