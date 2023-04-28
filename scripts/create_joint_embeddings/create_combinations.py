# import packages
from itertools import combinations
import pandas as pd
from scipy import stats

def per_feature_zscore(x):
    mean = x.mean()
    std = x.std(ddof=0)
    z = (x - mean)/std
    return z

embedding_path_list = snakemake.input['embedding_paths']
out_path = snakemake.output['out_path']

print("Normalization setting: " + snakemake.params['norm'])

embedding_df = []
for embedding_path in embedding_path_list:
    # load embeddings and add to list
    tmp_df = pd.read_csv(embedding_path, sep='\t')
    tmp_df = tmp_df.drop_duplicates(subset='gene_id').set_index('gene_id')
    #tmp_df = tmp_df.set_index('gene_id').transform(per_feature_zscore, axis=0)
    if snakemake.params['norm'] == "T":
        print('applying normalization')
        # zscore normalize embedding
        tmp_df = pd.DataFrame(stats.zscore(tmp_df.values, axis=None), 
                              index=tmp_df.index,
                              columns = tmp_df.columns)

    embedding_df.append(tmp_df)

# combine embeddings using the gene_id
embedding_df = pd.concat(embedding_df, join='inner', axis=1)
# set the gene_id as a column again
embedding_df = embedding_df.reset_index(drop=False)

# save to the given path
embedding_df.to_csv(out_path, sep='\t', index=False)
