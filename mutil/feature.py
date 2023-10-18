def get_feat_d_from_ID(feat_ID, mut_row):
    feat_d = dict()
    for d in mut_row["genomic features"]:
        if d["RegulonDB ID"] == feat_ID:
            feat_d = d
            break
    return feat_d


def get_feat_d(json, RegulonDB_ID=None, name=None):
    feat_d = dict()
    for d in json:
        if (RegulonDB_ID and d["RegulonDB ID"] == RegulonDB_ID) or (name and d["name"] == name):
            feat_d = d
            break
    return feat_d



# The below is currently only appropriate for genomic and genetic features 
# since it's not accounting for annotation types that don't have any entries
# for mutations.
import pandas as pd

def get_feat_mut_cnt_df_from_links(mut_df, feat_col_name, link_col_name):
    f_cnt_df = pd.DataFrame(columns=["length", "observed mutation count", "name"])
    for i, r in mut_df.iterrows():
        for f, links in r[link_col_name].items():
            f_d = get_feat_d(RegulonDB_ID=f, json=r[feat_col_name])
            if f_d["RegulonDB ID"] in f_cnt_df.index:
                f_cnt_df.loc[f_d["RegulonDB ID"], "observed mutation count"] += len(links)
            else:
                f_len = f_d["range"][1] - f_d["range"][0] + 1 
                df = pd.DataFrame({"name": f_d["name"], "length": f_len, "observed mutation count": len(links)},
                                    index=[f_d["RegulonDB ID"]])  # "name" column just for visual inspection
                f_cnt_df = f_cnt_df.append(df, sort=False)

    return f_cnt_df
