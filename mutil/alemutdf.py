import os
import os.path
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

module_path = os.path.abspath(os.path.join('.'))
if module_path not in sys.path:
    sys.path.append(module_path)


# TODO: Integrate functionality to get the series of fixed mutations defined by a final mutation with a frequency larger than a given threshold into this list of scripts. This functionality is implemented in substitution_rate_decline.ipynb.


def get_exp_afir(mut):
    return mut.exp + " " + get_afir(mut)


def get_afir(mut):
    return str(int(mut.ale)) + " " + str(int(mut.flask)) + " " + str(int(mut.isolate)) + " " + str(int(mut.tech_rep))


def get_coding_muts_df(mut_df):
    coding_mut_df = mut_df[~mut_df["Details"].str.contains("intergenic")
                           & ~mut_df["Details"].str.contains("pseudogene")
                           & ~mut_df["Details"].str.contains("noncoding")].copy()
    return coding_mut_df


def get_noncoding_muts_df(mut_df):
    coding_mut_df = mut_df[mut_df["Details"].str.contains("intergenic")
                           | mut_df["Details"].str.contains("pseudogene")
                           | mut_df["Details"].str.contains("noncoding")].copy()
    return coding_mut_df


def get_genetic_muts_df(mut_df):
    df = mut_df[~mut_df["Details"].str.contains("intergenic")].copy()
    return df


import traceback
import chardet


def get_mut_dataframe(CSV_file_path,
                      include_dups = False,
                      intragenic_muts_only = False):

    with open(CSV_file_path, 'rb') as rawdata:
        result = chardet.detect(rawdata.read(100000))

    # Step 1: Import database
    raw_db = pd.read_csv(CSV_file_path, encoding=result['encoding'])
    if 'Function' in raw_db.columns:
        raw_db = raw_db.drop('Function', axis=1)
    if 'Product' in raw_db.columns:
        raw_db = raw_db.drop('Product', axis=1)
    if 'GO Process' in raw_db.columns:
        raw_db = raw_db.drop('GO Process', axis=1)
    if 'GO Component' in raw_db.columns:
        raw_db = raw_db.drop('GO Component', axis=1)

    # TODO: The below will crash if most than these columns. Needs to ignore other mutation columns.
    # Step 2: Separate columns based on usage
    keep_cols = ['Position','Mutation Type','Sequence Change','Details','Gene']
    if "Reference Seq" in raw_db.columns:  # For backwards compatibility with older exported mutation data
        keep_cols.append("Reference Seq")
    if "Mut ID" in raw_db.columns:  # For backwards compatibility with older exported mutation data
        keep_cols.append("Mut ID")
    if 'Gene (Scrollable)' in raw_db.columns:  # For new ensemble ALEdb
        raw_db.rename(columns={'Gene (Scrollable)': 'Gene'}, inplace=True)
    mut_cols = sorted(list(set(raw_db.columns) - set(keep_cols)))

    file_name = os.path.basename(CSV_file_path)
    exp_name_from_file = file_name.replace(".csv", '')

    # Step 3: Shift mutation column names into row identifiers
    csv_file_mutat_df = pd.DataFrame()

    try:
        for col in mut_cols:
            df = raw_db[raw_db[col].notnull()][keep_cols]
            exp_name = '_'.join(col.split(' ')[:-4])
            if exp_name == '':  # Will happen with newer mutation data exported from ALEdb
                exp_name = exp_name_from_file
            df['exp'] = exp_name
            df['ale'] = int(col.split(' ')[-4][1:])
            df['flask'] = int(col.split(' ')[-3][1:])
            df['isolate'] = int(col.split(' ')[-2][1:])
            df['tech_rep'] = int(col.split(' ')[-1][1:])

            df['presence'] = raw_db[raw_db[col].notnull()][col]
            csv_file_mutat_df = pd.concat([csv_file_mutat_df,df])

        csv_file_mutat_df = csv_file_mutat_df[['exp','ale','flask','isolate','tech_rep','presence'] + keep_cols]
        csv_file_mutat_df = csv_file_mutat_df.fillna('')

        # Remove mutation entries with empty gene since they will screw up mutat_df.groupby(['Gene', ...])
        csv_file_mutat_df = csv_file_mutat_df.loc[csv_file_mutat_df['Gene'] != '']

        # Remove weird characters between gene names in multiple gene annotation.
        csv_file_mutat_df['Gene'] = csv_file_mutat_df['Gene'].str.replace("  ", " ")

        if not include_dups:
            csv_file_mutat_df = csv_file_mutat_df.loc[csv_file_mutat_df['Details'] != 'Duplication']

        if intragenic_muts_only:
            csv_file_mutat_df = csv_file_mutat_df.loc[csv_file_mutat_df['Gene'].str.contains(',') == False]
    except Exception as e:
        print("Issue with file: " + CSV_file_path)
        print(e)
        print(traceback.print_exc())

    return csv_file_mutat_df


def get_all_sample_mut_df(dir_path,
                          include_dups = False,
                          intragenic_muts_only = False):
    mutat_df = pd.DataFrame()
    mutat_df_list = []
    for file_name in os.listdir(dir_path):
        file_path = dir_path+'/'+file_name
        # print(file_path)
        mutat_df_list.append(get_mut_dataframe(file_path, include_dups, intragenic_muts_only))

    mutat_df = pd.concat(mutat_df_list)
    return mutat_df


def _get_exp_ale_set(mut_df):
    exp_ale_df = mut_df.copy()
    exp_ale_df["exp ale"] = exp_ale_df["exp"] + ' ' + exp_ale_df["ale"].map(str)
    exp_ale_set = set(exp_ale_df['exp ale'].tolist())
    return exp_ale_set


def get_gene_mut_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    column_to_delete_list = ["presence",
                         "tech_rep",
                         "isolate",
                         "flask",
                         "Position",
                         "Mutation Type",
                         "Sequence Change",
                         "Details"]
    current_columns = list(mut_df.columns.values)
    for column_to_delete in column_to_delete_list:
        if column_to_delete in current_columns:
            del mut_df[column_to_delete]
    mut_df = mut_df.drop_duplicates()
    mut_df["exp ale"] = mut_df["exp"]+' '+mut_df["ale"].map(str)

    mut_df_column_name_set = _get_exp_ale_set(mut_df)

    # Get mut_mat_df indexes of Gene names.
    # Can't simply use the "Gene" column of trunc_mut_df since a gene may be mutated in
    # more than one ALE exp, and we want a set of unique gene names.
    mut_df_index_set = set(mut_df["Gene"].tolist())

    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for gene_name, all_gene_mut_df in mut_df.groupby("Gene"):
        for index, mut_df_row in all_gene_mut_df.iterrows():
            mut_mat_df.loc[gene_name, mut_df_row["exp ale"]] = 1

    return mut_mat_df


def get_gene_mut_count_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    mut_df["exp ale"] = mut_df["exp"]+' '+mut_df["ale"].map(str)

    mut_df_column_name_set = _get_exp_ale_set(mut_df)

    # Get mut_mat_df indexes of Gene names.
    # Can't simply use the "Gene" column of trunc_mut_df since a gene may be mutated in
    # more than one ALE exp, and we want a set of unique gene names.
    mut_df_index_set = set(mut_df["Gene"].tolist())

    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for gene_name, all_gene_mut_df in mut_df.groupby("Gene"):
        for index, mut_df_row in all_gene_mut_df.iterrows():
            mut_mat_df.loc[gene_name, mut_df_row["exp ale"]] += 1

    return mut_mat_df


def get_mut_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    column_to_delete_list = ["presence",
                         "tech_rep",
                         "isolate",
                         "flask",
                         "Position",
                         "Mutation Type",
                         "Details"]
    current_columns = list(mut_df.columns.values)
    for column_to_delete in column_to_delete_list:
        if column_to_delete in current_columns:
            del mut_df[column_to_delete]
    mut_df = mut_df.drop_duplicates()
    mut_df["exp ale"] = mut_df["exp"] + ' ' + mut_df["ale"].map(str)
    mut_df["gene seq change"] = mut_df["Gene"] + ' ' + mut_df["Sequence Change"]

    mut_df_column_name_set = _get_exp_ale_set(mut_df)  # Get unique set of exp+ALE#

    mut_df_index_set = set(mut_df["gene seq change"].tolist())  # Get unique set of gene+(seq change) names.
    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for mut, all_mut_df in mut_df.groupby("gene seq change"):
        for index, mut_df_row in all_mut_df.iterrows():
            mut_mat_df.loc[mut, mut_df_row["exp ale"]] = 1

    return mut_mat_df


def get_enrichment_muts(mut_df):
    trunc_mut_df = mut_df.copy()

    # If we are going to keep Duplications, though we want to remove the '[' and ']' from their gene annotations.
#     trunc_mut_df["Gene"] = trunc_mut_df["Gene"].map(lambda x: x.lstrip('[').rstrip(']'))

    # Removing duplications
    trunc_mut_df = trunc_mut_df[trunc_mut_df["Details"] != "Duplication"]

    # Removing unused columns
    del trunc_mut_df["tech_rep"]
    del trunc_mut_df["isolate"]
    # Could have the same mutation, but with a different presence due to differences between clonal and population reseq'ing
    del trunc_mut_df["presence"]
    trunc_mut_df = trunc_mut_df.drop_duplicates()

    enrichment_mut_df = pd.DataFrame()
    for gene_mut_groupby in trunc_mut_df.groupby(["exp", "Gene"]):
        gene_mut_df = gene_mut_groupby[1]
        mutation_count = gene_mut_df.shape[0]
        if mutation_count > 1:
            enrichment_mut_df = enrichment_mut_df.append(gene_mut_df)

    return enrichment_mut_df


def get_ALE_final_flask_df(ALE_df):
    final_flask_num = ALE_df["flask"].max()
    final_flask_df = ALE_df[ALE_df["flask"]==final_flask_num]
    return final_flask_df


def get_exp_final_flask_df(exp_df):
    exp_final_flask_mut_df = pd.DataFrame()
    for ale, ale_mut_df in exp_df.groupby("ale"):
        ALE_final_flask_df = get_ALE_final_flask_df(ale_mut_df)
        exp_final_flask_mut_df = exp_final_flask_mut_df.append(ALE_final_flask_df)
    return exp_final_flask_mut_df


def _get_max_freq_mut_df(mut_df):
    max_mut_freq = mut_df["presence"].max()
    max_mut_df = mut_df[mut_df["presence"] == max_mut_freq]
    max_mut_df = max_mut_df.sort_values("flask")
    earliest_flask_max_mut_df = max_mut_df[:1]
    return earliest_flask_max_mut_df  # only need to return first row of max_mut_df


def get_ALE_max_freq_mut_df(ALE_mut_df, endpoint_flask_only):
    if endpoint_flask_only:
        ALE_mut_df = get_ALE_final_flask_df(ALE_mut_df)
    dfs = []
    #  List of everything that defines a mutation besides instance values (freq, isolate#, etc.)
    mutation_descriptors_list = ["Position", "Mutation Type", "Sequence Change", "Details", "Gene"]
    for mut_group in ALE_mut_df.groupby(mutation_descriptors_list):
        mut_df = mut_group[1]
        if len(mut_df) > 0:
            dfs.append(_get_max_freq_mut_df(mut_df))
    return pd.concat(dfs)


def get_exp_max_freq_mut_df(exp_mut_df, endpoint_flask_only):
    dfs = []
    for ale_name, ALE_mut_df in exp_mut_df.groupby(["ale"]):
        if len(ALE_mut_df) > 0:
            dfs.append(get_ALE_max_freq_mut_df(ALE_mut_df, endpoint_flask_only))
    return pd.concat(dfs)


def get_multi_exp_max_freq_mut_df(mut_df, endpoint_flask_only):
    dfs = []
    for group_name, exp_ale_mut_df in mut_df.groupby(["exp", "ale"]):
        if len(exp_ale_mut_df) > 0:
            dfs.append(get_ALE_max_freq_mut_df(exp_ale_mut_df, endpoint_flask_only))
    return pd.concat(dfs)


def get_filtered_mut_df(fixed_mut_df, filter_mut_df):
    if "flask" in filter_mut_df.columns:
        filter_mut_df = filter_mut_df.drop("flask", axis=1)
    if "isolate" in filter_mut_df.columns:
        filter_mut_df = filter_mut_df.drop("isolate", axis=1)
    if "tech_rep" in filter_mut_df.columns:
        filter_mut_df = filter_mut_df.drop("tech_rep", axis=1)
    if "presence" in filter_mut_df.columns:
        filter_mut_df = filter_mut_df.drop("presence", axis=1)
    return pd.merge(fixed_mut_df,
                    filter_mut_df,
                    how="inner",
                    on=["exp", "ale", "Position", "Mutation Type", "Sequence Change", "Details", "Gene"])


# Filtering for fixed mutation SERIES that have a max freq larger than mut_freq_floor.
def get_filtered_fixed_mut_series_df(fixed_mut_df, mut_freq_floor):

    # TODO: Seems like building the filter_mut_df could be replaced with _get_max_freq_mut_df(...),
    # though need to make unit test to clarify.
    filter_mut_df_list = []
    for exp_name, exp_mut_df in fixed_mut_df.groupby(["exp"]):
        filter_mut_df_list.append(get_exp_max_freq_mut_df(exp_mut_df, endpoint_flask_only=True))
    filter_mut_df = pd.concat(filter_mut_df_list)

    filter_mut_df = filter_mut_df[filter_mut_df["presence"]>mut_freq_floor]
    return get_filtered_mut_df(fixed_mut_df, filter_mut_df)

'''
# !!! get_ALE_mut_type_frac_df AND get_exp_avg_ale_mut_type_frac_df ARE OBSOLETE.
# USE get_mut_type_avg_frac_across_class_df INSTEAD

def get_ALE_mut_type_frac_df(mut_df, df_column_name):
    mut_type_set = set(mut_df[df_column_name])
    mut_type_frac_df = pd.DataFrame()
    for exp_ale_tuple, exp_ale_max_mut_df in mut_df.groupby(["exp", "ale"]):
        for mut_type in mut_type_set:
            mut_type_count = len(exp_ale_max_mut_df[exp_ale_max_mut_df[df_column_name]==mut_type])
            mut_type_fraction = mut_type_count/len(exp_ale_max_mut_df)

            mut_type_frac_df = mut_type_frac_df.append(pd.DataFrame([[exp_ale_tuple[0],
                                                                        exp_ale_tuple[1],
                                                                        mut_type,
                                                                        mut_type_fraction]],
                                                              columns=["exp",
                                                                       "ale",
                                                                       df_column_name,
                                                                       "fraction"]))
    return mut_type_frac_df


def get_exp_avg_ale_mut_type_frac_df(mut_df, column_name):
    ale_mut_type_frac_df = get_ALE_mut_type_frac_df(mut_df, column_name)
    exp_avg_ale_mut_type_frac_df = pd.DataFrame()
    for exp_name_mut_type_tuple, mut_type_df in ale_mut_type_frac_df.groupby(["exp", column_name]):
        exp_name = exp_name_mut_type_tuple[0]
        mut_type = exp_name_mut_type_tuple[1]
        mut_type_ale_mean = np.mean(mut_type_df["fraction"])
        df = pd.DataFrame([[exp_name, mut_type, mut_type_ale_mean]],
                          columns=["experiment", column_name, "fraction"])
        exp_avg_ale_mut_type_frac_df = exp_avg_ale_mut_type_frac_df.append(df)
    return exp_avg_ale_mut_type_frac_df    
'''


def get_mut_type_frac_across_class_df(mut_df, exp_level_l, mut_type_class, all_mut_type_class_set):
    mut_type_frac_df = pd.DataFrame()
    for class_tup, exp_ale_max_mut_df in mut_df.groupby(exp_level_l):
        for mut_type in all_mut_type_class_set:
            mut_type_count = len(exp_ale_max_mut_df[exp_ale_max_mut_df[mut_type_class]==mut_type])
            mut_type_frac = mut_type_count/len(exp_ale_max_mut_df)

            if len(exp_level_l) > 1:
                l = list(class_tup)
            else:
                l = [class_tup]
            data_l = l + [mut_type, mut_type_frac]
            col_name_l = exp_level_l + [mut_type_class, "fraction"]

            mut_type_frac_df = mut_type_frac_df.append(pd.DataFrame([data_l], columns=col_name_l))
    return mut_type_frac_df


def get_mut_type_avg_frac_across_class_df(mut_df, exp_level_l, mut_type_class, all_mut_type_class_set):
    mut_type_frac_df = get_mut_type_frac_across_class_df(mut_df, exp_level_l, mut_type_class, all_mut_type_class_set)
    class_avg_mut_type_frac_df = pd.DataFrame()

    # !!! Assuming the first class is the final class to consider penetration across
    final_penetration_class = exp_level_l[0]
    for tup, mut_type_df in mut_type_frac_df.groupby([final_penetration_class, mut_type_class]):
        final_class_type = tup[0]
        mut_type = tup[1]
        mut_type_class_type_mean = np.mean(mut_type_df["fraction"])
        df = pd.DataFrame([[final_class_type, mut_type, mut_type_class_type_mean]],
                          columns=[final_penetration_class, mut_type_class, "fraction"])
        class_avg_mut_type_frac_df = class_avg_mut_type_frac_df.append(df)

    return class_avg_mut_type_frac_df


def get_standardized_annotations(muts_df):
    # Different experiments have different strings for position (some with commas, some without),
    # therefore going ahead and changing them all to integers
    muts_df.Position = muts_df.Position.apply(lambda x: int(str(x).replace(",", "")))
    muts_df.Position = muts_df.Position.astype(int)

    # This work is also currently duplicated NB4. Keep the implementation here and remove the code from NB4
    muts_df["Gene"] = muts_df["Gene"].apply(lambda a: "rph" if a == "[rph], [rph]" else a)
    muts_df["Gene"] = muts_df["Gene"].apply(lambda a: "rph" if a == "[rph],[rph]" else a)
    muts_df.head()

    return muts_df