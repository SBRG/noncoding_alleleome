# Not true unit tests, but serves current purpose.
# Meant to be executed as a script using command line Python.

import pandas as pd
from alemutdf import get_multi_exp_max_freq_mut_df, get_all_sample_mut_df, get_mut_dataframe, get_ALE_final_flask_df, \
    get_exp_max_freq_mut_df, get_ALE_max_freq_mut_df, get_filtered_mut_df, get_coding_muts_df

MUT_DF_COLUMN_LIST = ['exp', 'ale', 'flask', 'isolate', 'tech_rep', 'presence', 'Position', 'Mutation Type',
                      'Sequence Change', 'Details', 'Gene']


# get_coding_muts_df unit test
test_mut_df = pd.DataFrame([
    ['K12 EXP2', 10, 247, 1, 1, 1.0, '111', 'SNP', 'G→T', 'intergenic (+435/‑21)', 'dapB, carA'],
    ['K12 EXP2', 10, 247, 1, 1, 1.0, '222', 'SNP', 'A→T', 'test_seq_change', 'geneW'],
    ['K12 EXP2', 10, 247, 1, 1, 1.0, '222', 'SNP', 'G→A', 'pseudogene', 'geneC'],
    ['K12 EXP2', 10, 247, 1, 1, 1.0, '228', 'SNP', 'C→A', 'noncoding (63/120 nt)', 'rrfH'],
    ['K12 EXP2', 10, 247, 1, 1, 1.0, '444', 'DEL', 'Δ82bp', '', '[rph], [rph]']],
    columns=MUT_DF_COLUMN_LIST)
expected_mut_df = pd.DataFrame([
    ['K12 EXP2', 10, 247, 1, 1, 1.0, '222', 'SNP', 'A→T', 'test_seq_change', 'geneW'],
    ['K12 EXP2', 10, 247, 1, 1, 1.0, '444', 'DEL', 'Δ82bp', '', '[rph], [rph]']],
    columns=MUT_DF_COLUMN_LIST)
output_mut_df = get_coding_muts_df(test_mut_df)
output_mut_df.reset_index(drop=True, inplace=True)  # Specific for unit test logic in finding equality between returned and expected.
assert(output_mut_df.equals(expected_mut_df))


# Test get_exp_max_freq_mut_df(...)
fixed_mut_df = get_mut_dataframe(CSV_file_path="test_data/test3.csv")
exp_max_freq_mut_df = get_exp_max_freq_mut_df(fixed_mut_df, endpoint_flask_only=True)
max_freq_mut_set = set(list(exp_max_freq_mut_df["presence"]))
assert(set([1, 0.9, 0.97, 0.95, 0.4]) == max_freq_mut_set)

# Test get_ALE_max_freq_mut_df(...)
# !!! "test_data/test2.csv" only has one ALE !!!
ALE_mut_df = get_mut_dataframe(CSV_file_path="test_data/test2.csv")
ALE_max_freq_mut_df = get_ALE_max_freq_mut_df(ALE_mut_df, endpoint_flask_only=False)
max_freq_mut_list = list(ALE_max_freq_mut_df["presence"])
assert(1 in max_freq_mut_list)
assert(0.9 in max_freq_mut_list)

ALE_max_freq_mut_df = get_ALE_max_freq_mut_df(ALE_mut_df, endpoint_flask_only=True)
max_freq_mut_list = list(ALE_max_freq_mut_df["presence"])
assert(1 in max_freq_mut_list)
assert(0.8 in max_freq_mut_list)


# Test for get_ALE_final_flask_df(...)
fixed_mut_df = get_mut_dataframe(CSV_file_path="test_data/test2.csv")
final_flask_df = get_ALE_final_flask_df(fixed_mut_df)
final_flask_num = final_flask_df.iloc[0]["flask"]
assert(final_flask_num==244)


# Test for get_max_freq_mut_df(...)
mut_df = get_all_sample_mut_df("test_data/exp_test_set/")
max_freq_mut_df = get_multi_exp_max_freq_mut_df(mut_df, endpoint_flask_only=True)
max_freq_mut_df = max_freq_mut_df.reset_index(drop=True)
expected_df = pd.DataFrame([['EXP1', 3, 244, 0, 2, 1.0, 1, 'SNP', 'C→A', 'P1100Q\xa0(CCG→CAG)\xa0', 'rpoB'],
 ['EXP1', 3, 244, 0, 2, 0.9, 2, 'SNP', 'T→G', 'intergenic\xa0', 'asdf, qwer'],
 ['EXP2', 3, 244, 0, 2, 1.0, 1, 'SNP', 'C→A', 'P1100Q\xa0(CCG→CAG)\xa0', 'rpoB'],
 ['EXP2', 3, 244, 0, 2, 0.9, 2, 'SNP', 'T→G', 'intergenic\xa0', 'asdf, qwer'],
 ['EXP2', 4, 3, 0, 1, 0.95, 3, 'DEL', 'Delta 3', '', 'qwer'],
 ['EXP2', 4, 3, 0, 1, 0.97, 4, 'SNP', 'C→T', 'D100Q\xa0(CCC→CCA)\xa0', 'fghj'],
 ['EXP2', 4, 3, 0, 1, 0.4, 5, 'SNP', 'A→G', 'B333G\xa0(TGG→GGG)\xa0', 'xcvb']], 
 columns=["exp","ale","flask","isolate","tech_rep","presence","Position","Mutation Type","Sequence Change","Details","Gene"])
assert(max_freq_mut_df.equals(expected_df) == True)


# Test for get_filtered_mut_df(...)
test_mut_df = get_mut_dataframe("test_data/get_filtered_fixed_mut_series_df.csv")
# simulating externally generated max freq filtered mutations 
filter_mut_df = get_mut_dataframe("test_data/get_filtered_fixed_mut_series_df.csv")
filter_mut_df = filter_mut_df.drop("flask", axis=1)
filter_mut_df = filter_mut_df.drop("isolate", axis=1)
filter_mut_df = filter_mut_df.drop("tech_rep", axis=1)
filter_mut_df = filter_mut_df.drop("presence", axis=1)
filter_mut_df = filter_mut_df.drop_duplicates()
filter_mut_df = filter_mut_df[filter_mut_df["Gene"] != "rpoB"]

filtered_fixed_mut_df = get_filtered_mut_df(test_mut_df, filter_mut_df)
expected_df = pd.DataFrame(
    [['C13', 3, 133, 1, 1, 1.0, 3815859, 'DEL', 'Δ82\xa0bp', '', '[rph], [rph]'],
     ['C13', 3, 57, 1, 1, 1.0, 3815859, 'DEL', 'Δ82\xa0bp', '', '[rph], [rph]'],
     ['C13', 3, 133, 1, 1, 1.0, 3815810, 'DEL', 'Δ1\xa0bp', 'intergenic\xa0(‑42/+24)', 'pyrE, rph'],
     ['C13', 3, 57, 1, 1, 1.0, 3815810, 'DEL', 'Δ1\xa0bp', 'intergenic\xa0(‑42/+24)', 'pyrE, rph'],
     ['C13', 5, 135, 1, 1, 1.0, 3815810, 'DEL', 'Δ1\xa0bp', 'intergenic\xa0(‑42/+24)', 'pyrE, rph'],
     ['C13', 5, 54, 1, 1, 1.0, 3815810, 'DEL', 'Δ1\xa0bp', 'intergenic\xa0(‑42/+24)', 'pyrE, rph']],
    columns=['exp', 'ale', 'flask', 'isolate', 'tech_rep', 'presence', 'Position', 'Mutation Type', 'Sequence Change', 'Details', 'Gene'])
assert(filtered_fixed_mut_df.equals(expected_df))


print("DONE")
