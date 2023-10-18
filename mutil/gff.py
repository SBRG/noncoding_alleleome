import pandas as pd


def group_sequences(sorted_l): 
    res = [[sorted_l[0]]] 
  
    for i in range(1, len(sorted_l)): 
        if sorted_l[i-1]+1 == sorted_l[i]: 
            res[-1].append(sorted_l[i]) 
  
        else: 
            res.append([sorted_l[i]]) 
    return res


def get_df_from_csv(csv_file):
    gff_df = pd.read_csv(
        csv_file,
        comment='#',
        header=None,
        sep='\t')
    d = {
        0:"seqname",
        1:"source",
        2:"feature",
        3:"start",
        4:"end",
        5:"score",
        6:"strand",
        7:"frame",
        8:"attribute"}
    gff_df = gff_df.rename(columns=d)
    return gff_df


def update_gff_feature_description(feature, attribute):
    updated_feature = feature
    if "Note=Substrate binding" in attribute:
        updated_feature = "Substrate binding site"
    if "Note=Allosteric FBP inhibitor binding" in attribute:
        updated_feature = "Allosteric FBP inhibitor binding site"
    if "Note=Zinc%3B shared with EIIA-Glc" in attribute:
        updated_feature = "EIIA and Zinc binding site"
    if "Note=Substrate;Ontology_term" in attribute:
        updated_feature = "Substrate binding site"
    if attribute == "Note=ATP":
        updated_feature = "ATP binding site"
    if "Note=ATP%3B via amide nitrogen" in attribute:
        updated_feature = "ATP binding site"
    if "Note=ATP%3B via carbonyl oxygen" in attribute:
        updated_feature = "ATP binding site"
    if "Note=HD;" in attribute:
        updated_feature = "HD domain"
    if "Note=TGS;" in attribute:
        updated_feature = "TGS domain"
    if "Note=ACT;" in attribute:
        updated_feature = "ACT domain"

    return updated_feature
