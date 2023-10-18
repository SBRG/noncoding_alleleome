# third-party packages
from Bio.SeqFeature import FeatureLocation
import pandas as pd

# =========================== FILE WRITING/READING ===========================
def write_fasta_lines(seq_tups, fpath):
    """
    Given a list of 2-tuples (seq_id, seq), write to a FASTA file at fpath
    """
    # make sure IDs are unique
    assert len(seq_tups) == len(set([tup[0] for tup in seq_tups]))
    with open(fpath, 'w') as f:
        f.write('\n'.join([f'>{seq_id}\n{insert_line_breaks(seq)}' for seq_id, seq in seq_tups]))


def insert_line_breaks(s, n_char_per_line=80):
    subs = []
    for i_line in range(0, len(s), n_char_per_line):
        subs.append(s[i_line:i_line+n_char_per_line])
    return '\n'.join(subs)


GFF_COL_NAMES = [
    'sequence_id', 'source', 'feature_type', 'left', 'right',
    'score', 'strand', 'phase', 'attributes'
]

def load_gff(gff_path):
    
    with open(gff_path, 'r') as f:
        current_line = f.readline()
        n_rows_to_skip = 0
        # skip the header
        while current_line[0] == '#' or current_line.strip() == '\n':
            n_rows_to_skip += 1
            current_line = f.readline()
        
        # now we can use pandas directly
        df = pd.read_csv(gff_path, skiprows=n_rows_to_skip, names=GFF_COL_NAMES, sep='\t')
        
    return df


# ============================= SEQUENCES =============================

# BioPython's handy FeatureLocation object ONLY accepts Python integers; cast to these
# make sure to handle the fact that bitome/GFFs uses 1-indexing
def get_sequence(extract_seq, left, right, strand):
    # handle the "wraparound" case - this would be a case where left is negative, or 
    # right is greater than seq length - the latter can easily be made the former
    # the range has circled all the way back around
    
    seq_len = len(extract_seq)
    
    # case where it's actually just all at the right end
    if left < 1 and right < 1:
        return FeatureLocation(int(seq_len + left - 1), int(seq_len + right), int(strand)).extract(extract_seq)
    # true wraparound cases
    # no way both are true - that would be more than the whole sequence
    elif left < 1 or right > seq_len:
        # these cases can be interconverted - consolidate
        if right > seq_len:
            left, right = left - seq_len, right - seq_len
        far_right_wraparound = FeatureLocation(
            int(seq_len + left - 1), seq_len, int(strand)
        ).extract(extract_seq)
        far_left_continued = FeatureLocation(
            0, int(right), int(strand)
        ).extract(extract_seq)
        # have to stitch the sequences together differently depending on the strand
        if strand == 1:
            return far_right_wraparound + far_left_continued
        else:
            return far_left_continued + far_right_wraparound
    # normal case
    else:
        return FeatureLocation(int(left)-1, int(right), int(strand)).extract(extract_seq)