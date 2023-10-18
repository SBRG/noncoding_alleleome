import pandas as pd


RNAP_GENE_NAMES = [
    "fecI",
    "fliA",
    "rpoH",
    "rpoA",
    "rpoB",
    "rpoC",
    "rpoN",
    "rpoD",
    "rpoE",
    "rpoS",
    "rpoF",
    "rpoZ"]


def get_coding_genetic_target_len_d(component_name_str, genes_df):
    genetic_target_len_d = dict()

    # Multiple gene mutation
    if ',' in component_name_str:
        # Get genes
        gene_l = component_name_str.split(',')
        for genetic_target in gene_l:
            if genetic_target in NON_REGULONDB_GENE_L:
                component_len = get_non_regulonDB_gene_len(genetic_target)
            else:
                component_len = len(genes_df[genes_df["GENE_NAME"] == genetic_target]["GENE_SEQUENCE"].iloc[0])
            genetic_target_len_d[genetic_target] = component_len
        # Get intergenic region
        idx = 0
        while idx < len(gene_l):
            # building intergenic region annotation
            intergenic_region_str = gene_l[idx] + '/' + gene_l[idx+1]
            idx+=2
            d = get_intergenic_len_d(intergenic_region_str, genes_df)
            if d[intergenic_region_str] != 0:
                genetic_target_len_d.update(d)
    elif ';' in component_name_str:  # The case of pseudogenes
        pesudogene_l = component_name_str.split(';')
        for genetic_target in pesudogene_l:
            component_len = len(genes_df[genes_df["GENE_NAME"] == genetic_target]["GENE_SEQUENCE"].iloc[0])
            genetic_target_len_d[genetic_target] = component_len
    else:
        if component_name_str in NON_REGULONDB_GENE_L:
            component_len = get_non_regulonDB_gene_len(component_name_str)
        else:
            component_len = len(genes_df[genes_df["GENE_NAME"] == component_name_str]["GENE_SEQUENCE"].iloc[0])
        genetic_target_len_d[component_name_str] = component_len
        
    return genetic_target_len_d


def get_intergenic_len_d(component_name_str, genes_df):
    intergenic_len = 0
    gene_l = component_name_str.split('/')
    gene_1_pos_d = get_gene_pos_d(gene_l[0], genes_df)
    gene_2_pos_d = get_gene_pos_d(gene_l[1], genes_df)
    if gene_1_pos_d["right"] == gene_2_pos_d["left"] or gene_2_pos_d["right"] == gene_1_pos_d["left"]:  
        intergenic_len = 0  # when no intergenic region between 2 genes
    elif gene_l[0] == gene_l[1]:  # same gene for whatever reason; see with "gatC/gatC" in 42C exp.
        intergenic_len = len(genes_df[genes_df["GENE_NAME"] == gene_l[0]]["GENE_SEQUENCE"].iloc[0])
    elif gene_1_pos_d["right"] < gene_2_pos_d["left"]:
        intergenic_len = gene_2_pos_d["left"] - gene_1_pos_d["right"] - 1
    else:
        intergenic_len = gene_1_pos_d["left"] - gene_2_pos_d["right"] - 1
    return {component_name_str: intergenic_len}


NON_REGULONDB_GENE_L = ["ykfN", "insZ", "ydbA", "yehH", "yhiS", "yjgX", "ybeM", "ykgP"]


def get_non_regulonDB_gene_len(gene_name):
    gene_len = 0
    if gene_name == "ykfN":
        gene_len = 63  # https://biocyc.org/gene?orgid=ECOLI&id=G0-16715
    elif gene_name == "insZ":
        gene_len = 897  # https://ecocyc.org/gene?orgid=ECOLI&id=G6632
    elif gene_name == "ydbA":
        gene_len = 6063  # https://ecocyc.org/gene?orgid=ECOLI&id=G6632
    elif gene_name == "yehH":
        gene_len = 2860  # https://biocyc.org/gene?orgid=ECOLI&id=G8208
    elif gene_name == "yhiS":
        gene_len = 1224  # https://biocyc.org/gene?orgid=ECOLI&id=G8208
    elif gene_name == "yjgX":
        gene_len = 1199  # https://biocyc.org/gene?orgid=ECOLI&id=G7898
    elif gene_name == "ybeM":
        gene_len = 788  # https://biocyc.org/gene?orgid=ECOLI&id=G7898
    elif gene_name == "ykgP":
        gene_len = 90  # https://biocyc.org/gene?orgid=ECOLI&id=G7898
    return gene_len


def get_non_regulonDB_gene_pos(gene_name):
    gene_pos_t = (0,0)
    if gene_name == "ykfN":
        gene_pos_t = (263150, 263212)
    elif gene_name == "insZ":
        gene_pos_t = (1294426, 1295322)
    elif gene_name == "ydbA":
        gene_pos_t = (1465392, 1474013)
    elif gene_name == "yehH":
        gene_pos_t = (2197410, 2200269)
    elif gene_name == "yhiS":
        gene_pos_t = (3651291, 3653713)
    elif gene_name == "yjgX":
        gene_pos_t = (4499593, 4500791)
    elif gene_name == "ybeM":
        gene_pos_t = (65831, 658818)
    elif gene_name == "ykgP":
        gene_pos_t = (313716, 313805)
    return gene_pos_t
    
    
def get_gene_pos_d(gene_name, genes_df):
    gene_pos_d = dict()
    if gene_name in NON_REGULONDB_GENE_L:
        gene_pos_t = get_non_regulonDB_gene_pos(gene_name)
        gene_pos_d = {
            "left": gene_pos_t[0],
            "right": gene_pos_t[1]
        }
    else:
        gene_pos_d = {"left": genes_df[genes_df["GENE_NAME"]==gene_name]["GENE_POSLEFT"].iloc[0],
                      "right": genes_df[genes_df["GENE_NAME"]==gene_name]["GENE_POSRIGHT"].iloc[0]}
    return gene_pos_d


# Essentially just expecting all_gene_bnum_name_df to be the local ./data/gene_name_syn_df.csv file.
# TODO: find a better way to load this data rather than expecting it as a parameter.
def get_gene_bnum(gene_name, all_gene_bnum_name_df):
    bnum = ''
    # all_gene_bnum_name_df = pd.read_csv('./data/gene_name_syn_df.csv')  # Too slow to do with every call
    # all_gene_bnum_name_df = pd.read_csv(gene_name_bnum_df_csv)
    gene_bnum_name_df = all_gene_bnum_name_df[all_gene_bnum_name_df["GENE_NAME"] == gene_name]
    if len(gene_bnum_name_df) > 0:
        bnum = gene_bnum_name_df.iloc[0]["BNUM"]
    return bnum


def make_regulondb_gene_synonym_df():
    gene_df = pd.read_csv("./data/gene.txt", sep="\t", comment='#', header=None)
    gene_df.columns = [
        "GENE_ID",
        "GENE_NAME",
        "GENE_POSLEFT",
        "GENE_POSRIGHT",
        "GENE_STRAND",
        "GENE_SEQUENCE",
        "GC_CONTENT",
        "CRI_SCORE",
        "GENE_NOTE",
        "GENE_INTERNAL_COMMENT",
        "KEY_ID_ORG",
        "GENE_TYPE",
        "GENE_NOTE_WEB"
    ]

    gene_synonym_df = pd.read_csv(
        "./data/object_synonym.txt",
        sep="\t",
        comment='#',
        header=None,
        quoting=3
    )
    gene_synonym_df.columns = ["OBJECT_ID", "OBJECT_SYNONYM_NAME", "OS_INTERNAL_COMMENT", "KEY_ID_ORG"]

    gene_id_name_df = pd.DataFrame()
    for _, g in gene_df.iterrows():
        srs = pd.Series({'GENE_ID': g['GENE_ID'], 'GENE_NAME': g['GENE_NAME']})
        df = pd.DataFrame(srs).T

        gn_syn_df = gene_synonym_df[gene_synonym_df["OBJECT_ID"] == g["GENE_ID"]]
        gn_syn_df = gn_syn_df[['OBJECT_ID', 'OBJECT_SYNONYM_NAME']]
        gn_syn_df.rename(columns={'OBJECT_ID': 'GENE_ID', 'OBJECT_SYNONYM_NAME': 'GENE_NAME'}, inplace=True)

        gene_id_name_df = pd.concat([gene_id_name_df, df, gn_syn_df])
    gene_id_name_df = gene_id_name_df[gene_id_name_df.GENE_NAME.notna()]

    gene_name_syn_df = pd.DataFrame()
    for g, gdf in gene_id_name_df.groupby(['GENE_ID']):
        bnum_df = gdf[gdf['GENE_NAME'].str.contains('^b\d{4}')]
        if len(bnum_df) > 0:
            bnum = bnum_df['GENE_NAME'].iloc[0]
            df = gdf.copy()
            df.rename(columns={'GENE_ID': 'BNUM'}, inplace=True)
            df.BNUM = bnum
            df = df[df['GENE_NAME'] != bnum]  # removing unnecessary bnum row
            gene_name_syn_df = pd.concat([gene_name_syn_df, df])
    gene_name_syn_df.BNUM = gene_name_syn_df.BNUM.apply(lambda s: s.replace(' (obsolete)', ''))

    gene_name_syn_df.to_csv('./data/gene_name_syn_df.csv', index=False)
