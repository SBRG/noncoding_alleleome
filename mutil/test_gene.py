import pandas as pd
from gene import get_gene_bnum

gene_name_syn_df = pd.read_csv('./aledbmutil/data/gene_name_syn_df.csv')
assert(get_gene_bnum("rpoB", gene_name_syn_df) == 'b3987')

print("DONE")
