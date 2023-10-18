import re

# Could enhance code structure through strategy pattern: https://sourcemaking.com/design_patterns/strategy, https://sourcemaking.com/design_patterns/strategy/python/1


CODON_NUCLEOTIDE_COUNT = 3

# HYPERMUTATOR_GENES contains all synonyms for genes known to manifest hypermutator phenotypes after being mutated.
# This list of genes is from 10.1016/S0966-842X(98)01424-3 and synonyms were retrieved from ecocyc.
HYPERMUTATOR_GENES = {
    "mutD", "dnaQ", "mutS", "ant", "plm", "fdv",
    "mutL", "mutH", "mutR", "topB", "prv", 'uvrD',
    "uvr502", "srjC", "uvrE", "dar-2", "dda",
    "mutU", "pdeB", "rad", "recL", "mutM",
    "fpg", "mutY", "micA", "mutT", "nudA",
    "dnaE", "polC", "sdgC", 'polA', "resA",
    "mutA", "glyV", "mutC", "glyW", "ins",
    "dam", "miaA", "trpX", "sodA", "sodB",
    "oxyR"  # We've had this one mutate in ALE experiments and not cause hypermutation, so don't include unless also checking if ALE is an outlier with mutation counts
    "nth", "nei", "xthA", "xth", "nfo",
    "ung", "vsr", "ada", "ogt", "recA",
    "zab", "umuB", "tif", "lexB", "recH", "rnmB", "srf",
    "recG", "radC", "spoV", "ssb", "exrB", "lexC"
                                           "hns"
    # We've had this one mutate in ALE experiments and not cause hypermutation, so don't include unless also checking if ALE is an outlier with mutation counts
}


def get_mutated_hypermutator_genes(m):
    mutated_hypermutator_genes = set()
    if m["coding"]:  # Mutation affect a gene's nucleotides
        mutated_genes_str = m["Gene"]
        mutated_genes = set(get_clean_mut_gene_list(mutated_genes_str))
        mutated_hypermutator_genes = mutated_genes & HYPERMUTATOR_GENES
    return mutated_hypermutator_genes


def get_MOB_type_str(MOB_seq_change_str):
    MOB_type_str = ""
    MOB_seq_change_str = MOB_seq_change_str.replace(u'\xa0', u' ')
    MOB_str_ID = ""
    if "IS" in MOB_seq_change_str:
        MOB_str_ID = "IS"
    elif "REP" in MOB_seq_change_str:
        MOB_str_ID = "REP"

    if MOB_str_ID != "":
        IS_str_start_idx = MOB_seq_change_str.find(MOB_str_ID)
        IS_str_end_idx = MOB_seq_change_str.find(' ', IS_str_start_idx)
        MOB_type_str = MOB_seq_change_str[IS_str_start_idx:IS_str_end_idx]

    return MOB_type_str


# If the DEL or INS is a multiple of 3, no matter what the position of the mutation within a frame,
# those frames downstream of the mutation will remaining in sync with the gene's coding frame.
# TODO: check if the remenants of combined frames result in a stop codon.
def is_frameshift(nuc_shift_size):
    return nuc_shift_size % 3 != 0


# TODO: pseudogenes are coding, therefore should update this
# and use a separate function for identifying pseudogenes.
def is_coding_mut(mut_details_str):
    is_coding = True
    noncoding_term_l = ["intergenic", "noncoding", "pseudogene"]
    if any(s in mut_details_str for s in noncoding_term_l):
        is_coding = False
    return is_coding


def is_genetic_mut(mut_details_str):
    is_genetic_mut = True
    nongenetic_term_l = ["intergenic"]
    if any(s in mut_details_str for s in nongenetic_term_l):
        is_genetic_mut = False
    return is_genetic_mut


def get_inv_size(seq_change_str):
    inv_size = 0
    seq_change_str = seq_change_str[:seq_change_str.find('bp')]
    inv_size = int(''.join([i for i in seq_change_str if i.isdigit()]))
    return inv_size


def get_con_size(seq_change_str):
    con_size = 0
    seq_change_str = seq_change_str[:seq_change_str.find('→')]
    con_size = int(''.join([i for i in seq_change_str if i.isdigit()]))
    return con_size


def get_MOB_size(seq_change_str):
    con_size = 0
    seq_change_str = seq_change_str[:seq_change_str.find('→')]
    con_size = int(''.join([i for i in seq_change_str if i.isdigit()]))
    return con_size


def get_sub_size(seq_change_str):
    sub_size = 0
    before_seq_change_str = seq_change_str[:seq_change_str.find('→')]
    before_seq_change_size = int(''.join([i for i in before_seq_change_str if i.isdigit()]))
    after_seq_change_str = seq_change_str[seq_change_str.find('→') + 1:]
    if "bp" in after_seq_change_str:
        after_seq_change_size = int(''.join([i for i in after_seq_change_str if i.isdigit()]))
        sub_size = before_seq_change_size - after_seq_change_size  # expecting these types of substitutions to always shrink in size
    else:  # no change in size, but multiple base pairs replaced.
        sub_size = before_seq_change_size
    return sub_size


def _get_before_after_seq_change_size(seq_change_str):
    before_after_seq_change_size = 0
    if '→' in seq_change_str:
        s = seq_change_str[:seq_change_str.find('→')]
        before_size_str = ''.join([i for i in s if i.isdigit()])
        s = seq_change_str[seq_change_str.find('→') + 1:]
        after_size_str = ''.join([i for i in s if i.isdigit()])
        before_after_seq_change_size = int(before_size_str) - int(after_size_str)
    return before_after_seq_change_size


def get_del_size(seq_change_str):
    del_size = 0
    if 'Δ' in seq_change_str or 'δ' in seq_change_str:
        del_size_str = ''.join([i for i in seq_change_str if i.isdigit()])
        del_size = int(del_size_str)
    if '→' in seq_change_str:
        del_size = _get_before_after_seq_change_size(seq_change_str)
    return del_size


def get_ins_size(seq_change_str):
    ins_size = 0
    if '→' in seq_change_str:
        before_seq_freq = int(seq_change_str[seq_change_str.find(')') + 1:seq_change_str.find('→')])
        after_seq_freq = int(seq_change_str[seq_change_str.find('→') + 1:])
        if "bp" in seq_change_str:
            seq_size = int(seq_change_str[seq_change_str.find('(') + 1:seq_change_str.find(' bp')])
        else:
            seq_str = seq_change_str[seq_change_str.find('(') + 1:seq_change_str.find(')')]
            seq_size = len(seq_str)
        ins_size = after_seq_freq * seq_size - before_seq_freq * seq_size
    if '+' in seq_change_str:
        if 'bp' in seq_change_str:
            ins_size = int(seq_change_str[seq_change_str.find('+') + 1:seq_change_str.find(' bp')])
        else:
            ins_size = len(seq_change_str[seq_change_str.find('+') + 1:])  # assuming it's just a sequence being added
    return ins_size


def get_amp_size(seq_change_str):
    seq_len_str = seq_change_str[:seq_change_str.find(' bp')]
    seq_len_str = seq_len_str.replace(',', '')
    seq_len = int(seq_len_str)
    multiplicity = int(seq_change_str[seq_change_str.find('x ') + len('x '):])
    seq_insertion_size = (multiplicity - 1) * seq_len
    return seq_insertion_size


def get_codon_pos_chng(codon_chng_str):
    ret_idx = 0
    codon_chng_list = codon_chng_str.split('→')
    for idx in range(CODON_NUCLEOTIDE_COUNT):  # assuming codon string is always going to be 3 long.
        if codon_chng_list[0][idx] != codon_chng_list[1][idx]:
            ret_idx = idx + 1
    return ret_idx


def is_in_TRN(details_str, gene_str, trn_gene_set):
    is_in_TRN = False
    if "intergenic" not in details_str:
        gene_list_str = gene_str
        gene_list_str = gene_list_str.replace('[', '')
        gene_list_str = gene_list_str.replace(']', '')
        gene_set = set(gene_list_str.split(', '))
        for gene in gene_set:
            if gene in trn_gene_set:
                is_in_TRN = True
    return is_in_TRN


# from https://en.wikipedia.org/wiki/Start_codon
def is_stop_codon(codon_str):
    codon_str = codon_str.lower()
    is_stop_codon = False
    stop_codon_l = ["taa", "tag", "tga"]
    if codon_str in stop_codon_l:
        is_stop_codon = True
    return is_stop_codon


# from https://en.wikipedia.org/wiki/Start_codon
def is_start_codon(codon_str):
    codon_str = codon_str.lower()
    is_start_codon = False
    start_codon_l = ["atg", "gtg", "ttg", "att", "ctg"]
    if codon_str in start_codon_l:
        is_start_codon = True
    return is_start_codon


def is_non_syn_SNP(amino_acid_change_str):
    return amino_acid_change_str[0] != amino_acid_change_str[-1]


def get_SNP_aa_pos(amino_acid_change_str):
    return int(re.sub("[^0-9]", "", amino_acid_change_str))


import math


def _replace_weirdly_encoded_dash_char(mut_details_str):
    weirdly_encoded_dash_char_1 = '‑'  # Breseq's weird encoding for the dash character
    weirdly_encoded_dash_char_2 = '‐'  # Another Breseq weirdly encoded dash character
    for c in [weirdly_encoded_dash_char_1, weirdly_encoded_dash_char_2]:
        if c in mut_details_str:
            mut_details_str = mut_details_str.replace(c, '-')        
    return mut_details_str


def get_DEL_INS_MOB_aa_start_pos(mut_details_str):
    aa_pos = None
    mut_details_str = _replace_weirdly_encoded_dash_char(mut_details_str)
    if len(mut_details_str):
        start_char = '('
        end_char = '/'
        if '-' in mut_details_str:
            end_char = '-'
        nuc_pos = int(mut_details_str[mut_details_str.find(start_char) + 1:mut_details_str.find(end_char)])
        aa_pos = math.ceil(nuc_pos / 3)
    return aa_pos


def get_DEL_INS_MOB_nuc_start_pos(mut_details_str):
    mut_details_str = _replace_weirdly_encoded_dash_char(mut_details_str)
    rel_nuc_pos = ''
    if len(mut_details_str):
        start_char = '('
        end_char = '/'
        if '-' in mut_details_str:
            end_char = '-'
        rel_nuc_pos = int(mut_details_str[mut_details_str.find(start_char) + 1:mut_details_str.find(end_char)])
    return rel_nuc_pos


def get_DEL_AA_range(mut_details_str):
    mut_details_str = _replace_weirdly_encoded_dash_char(mut_details_str)
    DEL_AA_range = ()
    if len(mut_details_str):

        start_char = '('
        end_char = '/'
        start_stop_str = mut_details_str[mut_details_str.find(start_char) + 1: mut_details_str.find(end_char)]
        l = start_stop_str.split('-')

        ints = [int(x) for x in l]
        if len(l) > 1:
            start_aa_pos = math.ceil(ints[0] / 3)
            stop_aa_pos = math.ceil(ints[1] / 3)
            DEL_AA_range = (start_aa_pos, stop_aa_pos)
        else:
            aa_pos_set = math.ceil(ints[0] / 3)
            DEL_AA_range = (aa_pos_set, aa_pos_set)
    return DEL_AA_range


def get_DEL_AA_set(mut_details_str):
    rng = get_DEL_AA_range(mut_details_str)
    return set(range(rng[0], rng[1] + 1))  # Have to add 1 since also want to consider final position in the range


def get_SUB_AA_range(mut_details_str):
    SUB_AA_range = ()
    mut_details_str = _replace_weirdly_encoded_dash_char(mut_details_str)
    if len(mut_details_str):
        start_char = '('
        end_char = '/'
        sub_str = mut_details_str[mut_details_str.find(start_char) + 1:mut_details_str.find(end_char)]
        SUB_AA_range = sub_str.split('-')
        SUB_AA_range[0] = int(SUB_AA_range[0])
        SUB_AA_range[1] = int(SUB_AA_range[1])
    return tuple(SUB_AA_range)


def get_codon_change_list(coding_SNP_details):
    codon_chng_str = get_codon_nuc_chng_str(coding_SNP_details)
    codon_change_list = codon_chng_str.split('→')
    return codon_change_list


def is_premature_stop_codon_SNP(coding_SNP_details):
    is_premature_stop_codon_SNP = False
    aa_chng_str = coding_SNP_details.split()[0]
    codon_chng_list = get_codon_change_list(coding_SNP_details)
    wt_codon = codon_chng_list[0]
    mut_codon = codon_chng_list[1]
    if is_non_syn_SNP(aa_chng_str) and is_stop_codon(mut_codon):
        is_premature_stop_codon_SNP = True
    return is_premature_stop_codon_SNP


def is_readthrough_codon_SNP(coding_SNP_details):
    is_readthrough_codon_SNP = False
    aa_chng_str = coding_SNP_details.split()[0]
    codon_change_list = get_codon_change_list(coding_SNP_details)
    wt_codon = codon_change_list[0]
    mut_codon = codon_change_list[1]
    if is_non_syn_SNP(aa_chng_str) \
            and is_stop_codon(wt_codon) and not is_stop_codon(mut_codon):
        is_readthrough_codon_SNP = True
    return is_readthrough_codon_SNP


def get_clean_mut_gene_list(gene_list_str):
    for s in STRINGS_TO_REMOVE:
        if s in gene_list_str:
            gene_list_str = gene_list_str.replace(s, "")
    if "genes" in gene_list_str:
        start_idx = gene_list_str.rfind("genes") + len("genes")
        gene_list_str = gene_list_str[start_idx:]

    split_str = ","
    if '|' in gene_list_str:
        split_str = "|"

    mut_gene_list = gene_list_str.split(split_str)
    clean_mut_gene_list = [gene for gene in mut_gene_list]
    return clean_mut_gene_list


def is_start_codon_removal(coding_SNP_details):
    is_start_codon_removal = False
    aa_chng_str = coding_SNP_details.split()[0]
    codon_change_list = get_codon_change_list(coding_SNP_details)
    wt_codon = codon_change_list[0]
    mut_codon = codon_change_list[1]
    if get_SNP_aa_pos(aa_chng_str) == 1 and is_start_codon(wt_codon) and not is_start_codon(mut_codon):
        is_start_codon_removal = True
    return is_start_codon_removal


STRINGS_TO_REMOVE = [' ', ">", "[", "]", "</i>", "<i>", "<b>", "</b>", "<BR>"]


def get_clean_mut_gene_list(gene_list_str):
    for s in STRINGS_TO_REMOVE:
        if s in gene_list_str:
            gene_list_str = gene_list_str.replace(s, "")
    if "genes" in gene_list_str:
        start_idx = gene_list_str.rfind("genes") + len("genes")
        gene_list_str = gene_list_str[start_idx:]

    split_str = ","
    if '|' in gene_list_str:
        split_str = "|"

    mut_gene_list = gene_list_str.split(split_str)
    clean_mut_gene_list = [gene for gene in mut_gene_list]
    return clean_mut_gene_list


def get_gene_count(mut_df_row):
    target_count = 0
    breseq_gene_annot_row = "Gene"
    if "mutation target annotation" in mut_df_row.keys():  # "mutation target annotation" used in AVA
        breseq_gene_annot_row = "mutation target annotation"
    gene_list_str = mut_df_row[breseq_gene_annot_row]
    mut_details_str = mut_df_row["Details"]
    if "intergenic" in mut_details_str:
        target_count = 1
    else:
        gene_set = set(get_clean_mut_gene_list(gene_list_str))
        target_count = len(gene_set)
    return target_count


STRUCTURAL_LEVEL = 0
OPERATIONAL_LEVEL = 1


# TODO: should check to see if start codons are ever destroyed.
def is_disruptive_SNP(mut_df_row):
    is_disruptive_SNP = False
    if mut_df_row["Mutation Type"].lower() == "snp" and mut_df_row["coding"]:
        if is_premature_stop_codon_SNP(mut_df_row["Details"]) or is_readthrough_codon_SNP(
                mut_df_row["Details"]) or is_start_codon_removal(mut_df_row["Details"]):
            is_disruptive_SNP = True
    return is_disruptive_SNP


# TODO: Not currently checking if pseudogenes are being further disrupted. Pseudogenes can still be translates; need to check if they get further disrupted.
def predict_mutation_effect_on_feature(mutation, feature):
    pred_eff = "other"

    if feature["feature type"] != "unknown":

        # Code block is for SVs
        if mutation["Mutation Type"].lower() in ["ins", "del", "mob"]:
            if feature["feature type"] == "gene":
                if is_frameshift(mutation["mutation size"]):
                    pred_eff = "truncation"
            elif mutation["mutation size"] >= 10:  # any feature (besides "unknown") if INS, DEL, or MOB > 10
                pred_eff = "truncation"

        # Code block is just for SNPs to genes
        if ((mutation["Mutation Type"].lower() == "snp") and (feature["feature type"] == "gene")):
            if is_disruptive_SNP(mutation):
                pred_eff = "truncation"
            else:
                aa_chng_str = mutation["Details"].split()[0]
                if is_non_syn_SNP(aa_chng_str):
                    pred_eff = "nonsynonymous"
                else:
                    pred_eff = "synonymous"

    return pred_eff


# For direct use with breseq mutations to a gene
def predict_mutation_effect_on_gene(mutation):
    pred_eff = "other"

    # Code block for SVs
    if mutation["Mutation Type"].lower() in ["ins", "del", "mob"]:
        if is_frameshift(mutation["mutation size"]):
            pred_eff = "truncation"
    
    # Code block for SNPs
    if mutation["Mutation Type"].lower() == "snp":
        if is_disruptive_SNP(mutation):
            pred_eff = "truncation"
        else:
            aa_chng_str = mutation["Details"].split()[0]
            if is_non_syn_SNP(aa_chng_str):
                pred_eff = "nonsynonymous"
            else:
                pred_eff = "synonymous"

    return pred_eff


def get_mob_size(seq_change_str, mob_sizes):
    mob_mut_size = 0

    s = seq_change_str
    s = s.replace("::", '')
    s = s.replace("(+)", '')
    s = s.replace("(–)", '')

    for annot_element in s.split("  "):
        if annot_element in mob_sizes.keys():
            mob_mut_size += mob_sizes[annot_element]
        else:
            s2 = annot_element
            s2 = s2.replace("bp", '')
            s2 = s2.replace(" ", '')
            s2 = s2.replace("Δ", '')
            s2 = s2.replace("+", '')

            val = 0
            if s2.isdigit():
                val = int(s2)
            else:
                val = len(s2)
            if 'Δ' in annot_element:
                val = (-1) * val

            mob_mut_size += val

    return mob_mut_size


# Currently not returning the size of MOBs. Isn't something currently necessary.
def get_mut_size(mut_df_row):
    mut_size = 0  # Currently defaulting everything except for INS and DEL to 0 since don't need them.
    if mut_df_row["Mutation Type"] == "SNP":
        mut_size = 1
    elif mut_df_row["Mutation Type"] == "INS":
        mut_size = get_ins_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "DEL":
        mut_size = get_del_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "INV":
        mut_size = get_inv_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "CON":
        mut_size = get_inv_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "SUB":
        mut_size = get_sub_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "AMP":
        mut_size = get_amp_size(mut_df_row["Sequence Change"])
    return mut_size


# Returns the range of mutations to nucleotides in mutation region before mutation.
def get_original_nuc_mut_range(mut_df_row):
    mut_range = (0, 0)
    if mut_df_row["Mutation Type"] == "SNP" \
            or mut_df_row["Mutation Type"] == "INS" \
            or mut_df_row["Mutation Type"] == "MOB" \
            or mut_df_row["Mutation Type"] == "AMP":
        mut_range = (mut_df_row["Position"], mut_df_row["Position"])
    elif mut_df_row["Mutation Type"] == "DEL" \
            or mut_df_row["Mutation Type"] == "INV" \
            or mut_df_row["Mutation Type"] == "CON" \
            or mut_df_row["Mutation Type"] == "SUB":
        mut_range = (mut_df_row["Position"], mut_df_row["Position"] - 1 + get_mut_size(mut_df_row))
    return mut_range


def get_codon_nuc_chng_str(coding_SNP_details):
    return coding_SNP_details[coding_SNP_details.find("(") + 1:coding_SNP_details.find(")")]


def get_coding_SNP_rel_nuc_pos(coding_SNP_details):
    rel_nuc_pos = ''
    codon_chng_str = get_codon_nuc_chng_str(coding_SNP_details)
    codon_nuc_pos = get_codon_pos_chng(codon_chng_str)
    aa_sub_str = coding_SNP_details[:coding_SNP_details.find(' ')]
    aa_sub = get_SNP_aa_pos(aa_sub_str)
    rel_nuc_pos = ((aa_sub - 1) * 3) + codon_nuc_pos
    return rel_nuc_pos


def get_genetic_coding_SNP_nuc_chng(SNP_details):
    nuc_chng = ''
    codon_change_list = get_codon_change_list(SNP_details)
    codon_chng_pos = get_codon_pos_chng(get_codon_nuc_chng_str(SNP_details))
    codon_chng_idx = codon_chng_pos - 1
    nuc_chng = codon_change_list[1][codon_chng_idx]
    return nuc_chng


def get_genetic_noncoding_or_pseudogene_SNP_nuc_chng(genetic_SNP_seq_change):
    return genetic_SNP_seq_change[-1]


def get_ins_seq(seq_change_str):
    ins_seq = ''
    if '→' in seq_change_str:
        before_seq_freq = int(seq_change_str[seq_change_str.find(')') + 1:seq_change_str.find('→')])
        after_seq_freq = int(seq_change_str[seq_change_str.find('→') + 1:])
        if "bp" in seq_change_str:  # TODO: find an example of this type of mutation; not currently sure what to expect with this on.
            # rare, more complicated to parse, and currently not occurring for mutation sets of interest, therefore not yet implementing.
            assert False, "needs to be implemented"
        else:
            single_ins_seq = seq_change_str[seq_change_str.find('(') + 1:seq_change_str.find(')')]
            ins_seq = single_ins_seq * (after_seq_freq - before_seq_freq)
    if '+' in seq_change_str:
        ins_seq = seq_change_str[seq_change_str.find('+') + 1:]
    return ins_seq
