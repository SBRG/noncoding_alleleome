from mut import get_inv_size, get_con_size, get_sub_size, get_del_size, get_ins_size, get_amp_size, is_coding_mut, \
    get_MOB_type_str, get_codon_pos_chng, is_premature_stop_codon_SNP, is_start_codon_removal, get_SNP_aa_pos, \
    get_gene_count, get_clean_mut_gene_list, get_DEL_INS_MOB_aa_start_pos, get_DEL_AA_set, \
    predict_mutation_effect_on_feature, get_mutated_hypermutator_genes, get_SUB_AA_range, get_DEL_AA_range, \
    get_DEL_INS_MOB_nuc_start_pos, get_coding_SNP_rel_nuc_pos, get_genetic_coding_SNP_nuc_chng, \
    get_genetic_noncoding_or_pseudogene_SNP_nuc_chng, get_ins_seq, is_readthrough_codon_SNP

assert (get_ins_seq("+ATT") == 'ATT')
assert (get_ins_seq('(T)6→7') == 'T')
assert (get_ins_seq('(T)6→9') == 'TTT')
assert (get_ins_seq('(AAAAGC)1→2') == 'AAAAGC')
assert (get_ins_seq('(AAAAGC)2→4') == 'AAAAGCAAAAGC')

assert (get_genetic_coding_SNP_nuc_chng('G157G (GGT→GGC)') == 'C')

assert (get_genetic_noncoding_or_pseudogene_SNP_nuc_chng('T→C') == 'C')

assert (get_coding_SNP_rel_nuc_pos(
    'A1B (AAA→TAA)') == 1)  # position 1 within codon (example codon given doesn't represent biological start codon)
assert (get_coding_SNP_rel_nuc_pos(
    'A1B (AAA→ATA)') == 2)  # position 2 within codon (example codon given doesn't represent biological start codon)
assert (get_coding_SNP_rel_nuc_pos(
    'A1B (AAA→AAT)') == 3)  # position 2 within codon (example codon given doesn't represent biological start codon)
assert (get_coding_SNP_rel_nuc_pos('A3B (AAA→TAA)') == 7)
assert (get_coding_SNP_rel_nuc_pos('A3B (AAA→ATA)') == 8)
assert (get_coding_SNP_rel_nuc_pos('A3B (AAA→AAT)') == 9)

assert (get_DEL_INS_MOB_nuc_start_pos('coding (798/1413 nt)') == 798)
assert (get_DEL_INS_MOB_nuc_start_pos('coding (50‐62/1203 nt)') == 50)

assert (get_mutated_hypermutator_genes({"coding": True, "Gene": "asdf, zxcv"}) == set())
assert (get_mutated_hypermutator_genes({"coding": False, "Gene": "asdf, zxcv, mutL"}) == set())
assert (get_mutated_hypermutator_genes({"coding": True, "Gene": "asdf, zxcv, mutL"}) == {"mutL"})
assert (get_mutated_hypermutator_genes({"coding": True, "Gene": "mutA, asdf, zxcv, mutL"}) == {"mutA", "mutL"})
assert (get_mutated_hypermutator_genes({"coding": True,
                                        "Gene": "asdf, zxcv, insB-5"}) == set())  # "ins" is a hypermutator genes; ensure that substrings aren't being identified as hypermutator genes with "insB-5".
assert (get_mutated_hypermutator_genes(
    {"coding": True, "Gene": '[ynaJ],uspE,fnr,ogt,abgT,abgB,abgA,abgR,mcaS,smrA,dgcM,zntB,fnrS,ynaL,dbpA,[ttcA]'}) == {
            "ogt"})  # testing complicated annotation string

i = "[araD],araA,araB"
e = ['araA', 'araB', 'araD']
assert (get_clean_mut_gene_list(i).sort() == e.sort())  # test without spaces between gene names
assert (get_clean_mut_gene_list('[rph],[rph]') == ['rph', 'rph'])
assert (get_clean_mut_gene_list("[asdf],zxcv,[qwer]") == ["asdf", "zxcv", "qwer"])
assert (get_clean_mut_gene_list("asdf, zxcv") == ["asdf", "zxcv"])  # test with space between gene names


def are_lists_equivalent(L1, L2):
    return len(L1) == len(L2) and sorted(L1) == sorted(L2)


gene_list_str = "geneW"
gene_list = ["geneW"]
output_list = get_clean_mut_gene_list(gene_list_str)
assert (are_lists_equivalent(output_list, gene_list))

gene_list_str = "[rph], [rph]"
gene_list = ["rph", "rph"]
output_list = get_clean_mut_gene_list(gene_list_str)
assert (are_lists_equivalent(output_list, gene_list))

gene_list_str = "sgrS|sgrT"
gene_list = ['sgrS', 'sgrT']
output_list = get_clean_mut_gene_list(gene_list_str)
assert (are_lists_equivalent(output_list, gene_list))

gene_list_str = "33 genesttcA>33 genes[ttcA], intR, ydaQ, ydaC, ralR, recT, recE, racC, ydaE, kilR, sieB, ydaF, ydaG, racR, ydaS, ydaT, ydaU, ydaV, ydaW, rzpR, rzoR, trkG, ynaK, ydaY, tmpR, lomR, insH1, lomR, stfR, tfaR, pinR, ynaE, > [ttcC]"
gene_list = ["ttcA", "intR", "ydaQ", "ydaC", "ralR", "recT", "recE", "racC", "ydaE", "kilR", "sieB", "ydaF", "ydaG",
             "racR", "ydaS", "ydaT", "ydaU", "ydaV", "ydaW", "rzpR", "rzoR", "trkG", "ynaK", "ydaY", "tmpR", "lomR",
             "insH1", "lomR", "stfR", "tfaR", "pinR", "ynaE", "ttcC"]
output_list = get_clean_mut_gene_list(gene_list_str)
assert (are_lists_equivalent(output_list, gene_list))

i = {"Gene": "[araD],araA,araB", "Details": ''}
assert (get_gene_count(i) == 3)

assert (is_coding_mut("intergenic (‑2/+1)") == False)
assert (is_coding_mut("A734V (GCG→GTG)") == True)
assert (is_coding_mut("noncoding (896/1542 nt)") == False)

assert (get_inv_size("1,443 bp inversion") == 1443)

assert (get_con_size("3024 bp→REL606:592057‑588495") == 3024)
assert (get_con_size("6 bp→REL606:1503176‑1503181") == 6)

assert (get_sub_size("57511 bp→6 bp") == 57505)
assert (get_sub_size('6827\xa0bp→81\xa0bp') == 6746)
assert (get_sub_size('2\xa0bp→CT') == 2)
assert (get_sub_size('2 bp→GA') == 2)

assert (get_del_size("Δ82 bp") == 82)
assert (get_del_size("δ82 bp") == 82)
assert (get_del_size("Δ1,199 bp") == 1199)
assert (get_del_size("(T)60→50") == 10)

assert (get_ins_size("(TTC)1→2") == 3)
assert (get_ins_size("+GCTA") == 4)
assert (get_ins_size("(45 bp)1→2") == 45)
assert (get_ins_size("+TG") == 2)
assert (get_ins_size("+46 bp") == 46)
assert (get_ins_size("+CGG") == 3)
assert (get_ins_size("(CGGTGGCTG)1→2") == 9)
assert (get_ins_size("(ATCG)1→4") == 12)

assert (get_amp_size('36,044 bp x 3') == 72088)

# assert(get_original_nuc_mut_range({"Mutation Type": "MOB", "Position": 99}) == (99,99))

assert (get_MOB_type_str("Δ1,199 bp") == "")
assert (get_MOB_type_str("IS1 (–) +8 bp") == "IS1")
assert (get_MOB_type_str("IS186 (–) +6 bp :: Δ1 bp") == "IS186")
assert (get_MOB_type_str("IS186 (–) +6 bp :: Δ1 bp") == "IS186")
assert (get_MOB_type_str("Δ1 :: IS186 (+) +6 bp :: Δ1") == "IS186")
test_str = "Δ6 bp :: REP161 (repetitive extragenic palindromic) element; contains 12 REP sequences (–) +2 bp :: +GGGGTGCCGCACTTCACAGCGGTGTAG"
assert (get_MOB_type_str(test_str) == "REP161")

assert (get_codon_pos_chng("ATC→AGC") == 2)
assert (get_codon_pos_chng("AAC→AAT") == 3)
assert (get_codon_pos_chng("AAC→AAC") == 0)

assert (is_premature_stop_codon_SNP('P1100Q\xa0(CCG→TAG)\xa0'))
assert (is_premature_stop_codon_SNP("R110G (CGT→GGT)") == False)
assert (is_premature_stop_codon_SNP("XXXAW (TAT→TAA)") == True)

assert (is_readthrough_codon_SNP("R110G (CGT→GGT)") == False)
assert (is_readthrough_codon_SNP("XXXAW (TAA→TAC)") == True)

assert (get_SNP_aa_pos("P1100Q") == 1100)

assert (is_start_codon_removal('P2Q\xa0(CTG→CCG)\xa0') == False)
assert (is_start_codon_removal('P1100Q\xa0(→TAG)\xa0') == False)
assert (is_start_codon_removal('M1M (ATG→ATA) †'))

# dictionaries can be substituted for Pandas series (DataFrame rows).
# testing pseudogene INS, DEL, MOB logic
m = {"Mutation Type": "DEL", "coding": False, "mutation size": 1}
f = {"genetic": True, "feature type": "pseudogene"}
assert (predict_mutation_effect_on_feature(m, f) == "other")

# start codon removal
m = {"Mutation Type": "SNP", "coding": True, "Details": "M1I (ATG→ATA) †"}
f = {"genetic": True, "feature type": "gene"}
assert (predict_mutation_effect_on_feature(m, f) == "truncation")

# premature stop codon
m = {"Mutation Type": "SNP", "coding": True, "Details": "P1100Q (CCG→TAG)"}
f = {"genetic": True, "feature type": "gene"}
assert (predict_mutation_effect_on_feature(m, f) == "truncation")

m = {"Mutation Type": "SNP", "coding": True, "Details": "A734V (GCG→GTG)"}
f = {"genetic": True, "feature type": "gene"}
assert (predict_mutation_effect_on_feature(m, f) == "nonsynonymous")

m = {"Mutation Type": "SNP", "coding": True, "Details": "A734A (GCG→GTG)"}
f = {"genetic": True, "feature type": "gene"}
assert (predict_mutation_effect_on_feature(m, f) == "synonymous")

for mt in ["AMP", "CNV", "SUB", "INV"]:
    assert (predict_mutation_effect_on_feature({"Mutation Type": "mt"}, {"feature type": "gene"}) == "other")

m = {"Mutation Type": "INS", "mutation size": 3}
f = {"genetic": True, "feature type": "gene"}
assert (predict_mutation_effect_on_feature(m, f) == "other")

m = {"Mutation Type": "INS", "mutation size": 2}
f = {"genetic": True, "feature type": "gene"}
assert (predict_mutation_effect_on_feature(m, f) == "truncation")

m = {"Mutation Type": "INS", "mutation size": 2}
f = {"genetic": False, "feature type": "promoter"}
assert (predict_mutation_effect_on_feature(m, f) == "other")

m = {"Mutation Type": "INS", "mutation size": 10}
f = {"genetic": False, "feature type": "TFBS"}
assert (predict_mutation_effect_on_feature(m, f) == "truncation")

# testing output for "unknown" (intergenic region)
f = {"feature type": "unknown"}
assert (predict_mutation_effect_on_feature(None, f) == "other")

weirdly_encoded_dash_char_1 = '‑'  # Breseq's weird encoding for the dash character
weirdly_encoded_dash_char_2 = '‐'  # Another Breseq weirdly encoded dash character
assert (get_DEL_INS_MOB_aa_start_pos("coding (1/20 nt)") == 1)
assert (get_DEL_INS_MOB_aa_start_pos("coding (4/20 nt)") == 2)
assert (get_DEL_INS_MOB_aa_start_pos("coding (58"+weirdly_encoded_dash_char_1+"61/1413 nt)") == 20)  # one type of weird hash character
assert (get_DEL_INS_MOB_aa_start_pos("coding (58"+weirdly_encoded_dash_char_2+"61/1413 nt)") == 20)  # another type of weird hash character

assert (get_DEL_AA_range('coding (50‐62/1203 nt)') == (17, 21))

assert (get_DEL_AA_set("coding (1/20 nt)") == {1})  # test the removal of 1 BP from AA1
assert (get_DEL_AA_set("coding (2/20 nt)") == {1})  # test the removal of 1 BP from AA1
assert (get_DEL_AA_set("coding (3/20 nt)") == {1})  # test the removal of 1 BP from AA1
assert (get_DEL_AA_set("coding (4/20 nt)") == {2})
assert (get_DEL_AA_set("coding (1‑2/20 nt)") == {1})  # test removal of BPs within AA1
assert (get_DEL_AA_set("coding (1‑2/20 nt)") == {1})  # test removal of BPs within AA1
assert (get_DEL_AA_set("coding (1‑3/20 nt)") == {1})  # test removal of BPs within AA1
assert (get_DEL_AA_set("coding (4‑6/20 nt)") == {2})  # test of removal of all BPs for one AA
assert (get_DEL_AA_set("coding (4‑9/20 nt)") == {2, 3})  # test for removal of all BPs for multiple AAs.
assert (get_DEL_AA_set("coding (6‑9/20 nt)") == {2, 3})  # test for removal of some BPs for multiple AAs.
assert (get_DEL_AA_set("coding (6‑8/20 nt)") == {2, 3})  # test for removal of some BPs for multiple AAs.
assert (get_DEL_AA_set("coding (4‑11/20 nt)") == {2, 3, 4})  # test for removal of some BPs for 3 AAs.
assert (get_DEL_AA_set("coding (2412‑2420/2547 nt") == {804, 805, 806,
                                                        807})  # 2412 is divisible by 3 == final nuc in an AA.

assert (get_SUB_AA_range('coding (114‐115/1260 nt)') == (114, 115))

print("DONE")
