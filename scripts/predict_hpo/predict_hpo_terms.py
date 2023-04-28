#!/usr/bin/env python
# coding: utf-8

import numpy as np
import polars as pl

import torch
from torch.utils.data import TensorDataset, DataLoader
import torch.nn.functional as F
import torch.nn as nn

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
import matplotlib.pyplot as plt


hpo = pl.read_csv(snakemake.input['gene2pheno'], sep = '\t')

# filter for hpo terms that have at least 20 genes annotated.
hpo = hpo.filter(pl.col('n_genes_annotated') >= 20)
hpo = hpo.with_columns(pl.lit(1).alias('one'))

hpo.filter(pl.col('HPO_Term_ID')== "HP:0012166")

hpo_labels = hpo.pivot(values="one", index="Ensembl_ID", columns="HPO_Term_ID")
hpo_labels = hpo_labels.fill_null(0)
hpo_labels

# read embedding
# embedding_path = '../data/AE_Embeddings/combinations/pca-gtex_pca-crispr_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall.tsv'
embedding_path = snakemake.input['emb']

emb = pl.read_csv(embedding_path, sep = '\t')
emb.shape

labels = emb.select('gene_id').join(hpo_labels, left_on='gene_id', right_on = 'Ensembl_ID', how = 'left').fill_null(0)


ds = TensorDataset(torch.as_tensor(emb.drop('gene_id').to_numpy().astype('float32')), 
                   torch.as_tensor(labels.drop('gene_id').to_numpy().astype('float32')))


# create train test indices.
np.random.seed(1234)
indices = np.arange(emb.shape[0])
np.random.shuffle(indices)
train_idx, test_idx = indices[:int(0.7 * emb.shape[0])], indices[int(0.7 * emb.shape[0]): ]

train_idx 
test_idx

train_ds = torch.utils.data.Subset(ds, train_idx)
test_ds = torch.utils.data.Subset(ds, test_idx)

train_dl = DataLoader(train_ds, batch_size=20_000)
test_dl = DataLoader(test_ds, batch_size=20_000)

x, y = next(iter(train_dl))
x_test, y_test = next(iter(test_dl))



class Model(torch.nn.Module):
    def __init__(self, input_size, hidden_size, out_size):
        super().__init__()
        
        #layers
        self.l1 = nn.Linear(input_size, hidden_size)
        self.l2 = nn.Linear(hidden_size, out_size)
    
    def forward(self, x):
        #one forward pass
        out = F.relu(self.l1(x))
        out = self.l2(out)
        return out


mod = Model(x.shape[1], 1000, y.shape[1])

loss_fn = nn.BCEWithLogitsLoss()
opt = torch.optim.Adam(mod.parameters())

for i in range(500):
    opt.zero_grad()
    loss = loss_fn(mod(x), y)
    loss.backward()
    opt.step()
    with torch.no_grad():
        test_loss = loss_fn(mod(x_test), y_test)
    print(loss, test_loss)

with torch.no_grad():
    y_pred_test = mod(x_test)


import pandas as pd
    
# write predictions to file for evaluation.
test_pred = pd.DataFrame(y_pred_test,
                         columns = labels.columns[1:], 
                         index = labels.select("gene_id")[test_idx].to_numpy().flatten()
                        )
test_pred.index.name = 'gene_id'
test_pred.to_csv(snakemake.output['hpo_pred'], sep = '\t')

    


# fig, axs = plt.subplots(4, 4, figsize=(24, 24))

# np.random.seed(59)
# for i, idx in enumerate(np.random.randint(y_pred_test.shape[1], size = 16)):
#     prec, recall, _ = precision_recall_curve(y_test.numpy()[:,idx], y_pred_test[:,idx])
#     pr_display = PrecisionRecallDisplay(precision=prec, recall=recall, pos_label=labels.columns[1:][idx] )
#     pr_display.plot(ax = axs[int(np.floor(i/4)), i%4])


# # In[36]:


# from sklearn.metrics import precision_recall_curve
# from sklearn.metrics import PrecisionRecallDisplay

# fig, axs = plt.subplots(1, figsize=(8, 8))

# prec, recall, _ = precision_recall_curve(y_test.numpy().flatten(), y_pred_test.flatten())
# pr_display = PrecisionRecallDisplay(precision=prec, recall=recall)
# pr_display.plot(ax = axs)


# # In[37]:


# from sklearn.metrics import average_precision_score


# # In[38]:


# average_precision_score_list = []
# for i in range(y_pred_test.shape[1]):
#     average_precision_score_list.append(
#         average_precision_score(y_test.numpy()[:,i], y_pred_test[:,i])
#     )


# # In[39]:


# len(average_precision_score_list)


# # In[40]:


# len(labels.columns[1:])


# # In[41]:


# num_genes_per_term = labels.sum().transpose(include_header=True).drop_nulls()
# num_genes_per_term.columns = ['hpo_term', 'n_genes']
# num_genes_per_term = num_genes_per_term.with_column(pl.col('n_genes').cast(pl.Int32))
# num_genes_per_term


# # In[42]:


# precision_scores = pl.DataFrame({"hpo_term": labels.columns[1:],
#               "average_precision_score": average_precision_score_list})
# precision_scores = precision_scores.sort('average_precision_score', reverse = True)
# precision_scores


# # In[43]:


# precision_scores = precision_scores.join(num_genes_per_term, on='hpo_term')


# # In[44]:


# import plotnine as p9


# # In[45]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes('hpo_term', 'average_precision_score')) + p9.geom_bar(stat = 'identity') + p9.theme(figure_size=(24, 7), axis_text_x=p9.element_text(rotation=90, hjust=1))


# # In[46]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes('n_genes', 'average_precision_score')) + p9.geom_point() + p9.scale_x_log10() + p9.theme_minimal() #+ p9.theme(figure_size=(24, 7), axis_text_x=p9.element_text(rotation=90, hjust=1))


# # In[47]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes('n_genes', 'average_precision_score')) + p9.geom_bin2d(bins = 50) + p9.scale_x_log10() + p9.theme_minimal() #+ p9.theme(figure_size=(24, 7), axis_text_x=p9.element_text(rotation=90, hjust=1))


# # In[ ]:





# # In[48]:


# # read HPO filler as a comparison


# # In[49]:


# HPOFiller = pl.read_csv('/data/ouga/home/ag_gagneur/brechtma/Downloads/prediction_HPOFiller_temporalpa_20190212_20200608.txt',
#                         sep = '\t', has_header = False, new_columns = ['rank', 'UniProt', 'HPO'])


# # In[50]:


# mapping = pl.read_csv('/data/ouga/home/ag_gagneur/brechtma/Downloads/mart_export_UNIProt_mapping.txt', sep = '\t')
# mapping = mapping.select([pl.col("Gene stable ID").alias("gene_id"), pl.col("UniProtKB/Swiss-Prot ID").alias("UniProt")])
# mapping = mapping.drop_nulls().unique()
# print(mapping.shape)

# HPOFiller = HPOFiller.join(mapping, on='UniProt')
# HPOFiller


# # In[51]:


# genes_in_test = labels.select("gene_id")[test_idx]


# # In[52]:


# hpo_rank = HPOFiller.pivot(values = "rank", index = "gene_id", columns = "HPO")


# # In[53]:


# hpo_rank = genes_in_test.join(hpo_rank, on = 'gene_id', how = "left")
# hpo_rank


# # In[54]:


# hpo_rank.sum()


# # In[55]:


# np.array(labels.columns)[np.isin(labels.columns, hpo_rank.columns)]


# # In[56]:


# labels.columns[5]


# # In[57]:


# fig, axs = plt.subplots(4, 4, figsize=(24, 24))

# np.random.seed(59)
# for i, idx in enumerate(np.random.randint(y_pred_test.shape[1], size = 16)):
    
#     # get Hpo name
#     hpo_term = labels.columns[1:][idx]
    
    
#     prec, recall, _ = precision_recall_curve(y_test.numpy()[:,idx], -hpo_rank.select(hpo_term).fill_null(1e9).to_numpy())
#     pr_display = PrecisionRecallDisplay(precision=prec, recall=recall, pos_label= hpo_term)
#     pr_display.plot(ax = axs[int(np.floor(i/4)), i%4], label = "HPO_Filler")


# # In[58]:


# hpo_rank.select(labels.columns[5]).fill_null(1e9).to_numpy()


# # In[59]:


# fig, axs = plt.subplots(4, 4, figsize=(24, 24))

# np.random.seed(59)
# for i, idx in enumerate(np.random.randint(y_pred_test.shape[1], size = 16)):
    
#      # get Hpo name
#     hpo_term = labels.columns[1:][idx]
    
    
#     prec, recall, _ = precision_recall_curve(y_test.numpy()[:,idx], -hpo_rank.select(hpo_term).fill_null(1e9).to_numpy())
#     pr_display = PrecisionRecallDisplay(precision=prec, recall=recall, pos_label= hpo_term)
#     pr_display.plot(ax = axs[int(np.floor(i/4)), i%4], label = "HPO_Filler")

    
#     prec, recall, _ = precision_recall_curve(y_test.numpy()[:,idx], y_pred_test[:,idx])
#     pr_display = PrecisionRecallDisplay(precision=prec, recall=recall, pos_label=labels.columns[1:][idx] )
#     pr_display.plot(ax = axs[int(np.floor(i/4)), i%4], label = "Embedding")
    
#     axs[int(np.floor(i/4)), i%4].legend(loc = "upper right")
    


# # In[104]:


# fig, axs = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
# #fig.wide_layout()

# size = 25
# SMALL_SIZE = size
# MEDIUM_SIZE = size
# BIGGER_SIZE = size

# plt.rcParams.update({'font.size': MEDIUM_SIZE})
# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


# np.random.seed(59)
# for i, idx in enumerate(np.random.randint(y_pred_test.shape[1], size = 16)[np.array([2,3,6,7])]):
    
#      # get Hpo name
#     hpo_term = labels.columns[1:][idx]
    
    
#     prec, recall, _ = precision_recall_curve(y_test.numpy()[:,idx], -hpo_rank.select(hpo_term).fill_null(1e9).to_numpy())
#     pr_display = PrecisionRecallDisplay(precision=prec, recall=recall, pos_label= None)
#     pr_display.plot(ax = axs[int(np.floor(i/2)), i%2], label = "HPO_Filler")

    
#     prec, recall, _ = precision_recall_curve(y_test.numpy()[:,idx], y_pred_test[:,idx])
#     pr_display = PrecisionRecallDisplay(precision=prec, recall=recall, pos_label=None )
#     pr_display.plot(ax = axs[int(np.floor(i/2)), i%2], label = "Embedding")
    
#     axs[int(np.floor(i/2)), i%2].legend(loc = "upper right")
#     axs[int(np.floor(i/2)), i%2].title.set_text(hpo_term)


# # In[61]:


# def compute_HPO_Filler_precision():
#     try:
#         return average_precision_score(y_test.numpy()[:,i], -hpo_rank.select(hpo_term).fill_null(1e9).to_numpy())
#     except:
#         return np.nan
    

# average_precision_score_list_emb = []
# average_precision_score_list_HPO_Filler = []
# for i in range(y_pred_test.shape[1]):
    
#     # get Hpo name
#     hpo_term = labels.columns[1:][i]
    
#     average_precision_score_list_emb.append(
#         average_precision_score(y_test.numpy()[:,i], y_pred_test[:,i])
#     )
    
#     average_precision_score_list_HPO_Filler.append(
#         compute_HPO_Filler_precision()
#     )
    
    


# # In[62]:


# len(average_precision_score_list_HPO_Filler)


# # In[63]:


# len(labels.columns[1:])


# # In[64]:


# len(average_precision_score_list_emb)


# # In[65]:


# precision_scores = pl.DataFrame({"hpo_term": labels.columns[1:],
#               "average_precision_score_embedding": average_precision_score_list_emb, 
#               "average_precision_score_list_HPO_Filler": average_precision_score_list_HPO_Filler})
# precision_scores = precision_scores.sort('average_precision_score_embedding', reverse = True)
# precision_scores


# # In[66]:


# precision_scores = precision_scores.join(num_genes_per_term, on='hpo_term')


# # In[67]:


# precision_scores.mean()


# # In[68]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes('average_precision_score_list_HPO_Filler', 'average_precision_score_embedding')) + p9.geom_point(size = 0.5) + p9.theme_minimal() + p9.scale_x_log10() + p9.scale_y_log10() + p9.geom_abline()#+ p9.theme(figure_size=(24, 7), axis_text_x=p9.element_text(rotation=90, hjust=1))


# # In[69]:


# precision_scores = precision_scores.with_column(pl.when(pl.col("n_genes") < 50).then(pl.lit("< 50")).otherwise(pl.lit("> 50")).alias("bin"))


# # In[70]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes('average_precision_score_list_HPO_Filler', 'average_precision_score_embedding')) + p9.geom_point(size = 0.5) + p9.theme_minimal() + p9.scale_x_log10() + p9.scale_y_log10() + p9.geom_abline() + p9.facet_wrap("bin")


# # In[71]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes('n_genes', 'average_precision_score_embedding')) + p9.geom_point() + p9.scale_x_log10() + p9.theme_minimal() #+ p9.theme(figure_size=(24, 7), axis_text_x=p9.element_text(rotation=90, hjust=1))


# # In[72]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes('n_genes', 'average_precision_score_list_HPO_Filler')) + p9.geom_point() + p9.scale_x_log10() + p9.theme_minimal() #+ p9.theme(figure_size=(24, 7), axis_text_x=p9.element_text(rotation=90, hjust=1))


# # In[73]:


# # check some not so great examples


# # In[74]:


# precision_scores.filter(
#     (pl.col("n_genes") > 30) & 
#     (pl.col("n_genes") < 50) &
#     (pl.col("average_precision_score_embedding") < 0.01) & 
#     (pl.col("average_precision_score_list_HPO_Filler") > 0.1)
# )


# # In[75]:


# labels.select(["gene_id", "HP:0001765"]).filter(pl.col("HP:0001765") == 1)


# # In[76]:


# labels.select("HP:0001765" )[train_idx].sum()


# # In[77]:


# labels.select("HP:0001765" )[test_idx].sum()


# # In[78]:


# precision_scores.filter(pl.col("hpo_term") == "HP:0001780")


# # In[79]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes("n_genes")) + p9.geom_histogram(bins = 100) + p9.theme_minimal() + p9.scale_x_log10()


# # In[80]:


# precision_scores.describe()


# # In[81]:


# import pandas as pd


# # In[82]:


# precision_scores['bin'] = pd.qcut(precision_scores.n_genes, q=5).astype('str')


# # In[83]:


# precision_scores


# # In[84]:


# precision_scores.groupby('bin').mean()


# # In[85]:


# precision_scores.with_columns([pl.col('average_precision_score_list_HPO_Filler').mean().over('bin').alias('mean_HPO_Filler'),
#                               pl.col('average_precision_score_list_HPO_Filler').mean().over('bin').alias('mean_HPO_Filler')])


# # In[86]:


# p9.ggplot(precision_scores.to_pandas(), p9.aes('average_precision_score_list_HPO_Filler', 'average_precision_score_embedding')) + \
#     p9.geom_point(size = 0.5) + \
#     p9.geom_point(data = precision_scores.groupby('bin').mean().to_pandas(), color = 'red') + \
#     p9.theme_minimal() + \
#     p9.scale_x_log10() + p9.scale_y_log10() + p9.geom_abline() + p9.facet_wrap("bin", nrow = 1) + \
#     p9.theme(figure_size=(15, 8))


# # In[ ]:




