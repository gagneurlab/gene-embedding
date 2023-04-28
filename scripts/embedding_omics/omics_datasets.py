


import torch
import numpy as np
import polars as pl
import pandas as pd
import umap

import scipy.stats

from sklearn.decomposition import PCA
from torch.utils.data import Dataset, DataLoader, BatchSampler, RandomSampler, random_split
import plotnine as p9
import anndata
from scipy.sparse import coo_array

# insert path to project here:
project_path = '...path_to_project...'

  
class FasterDepMapDataset(Dataset):
    def __init__(self,
                 path = project_path + 'dep_map_data/DepMap18Q3_gene_effect_mapped.tsv',
                 gene_index_path = project_path + '/processed_data/omics_index/gene_index.tsv'):
        
        self.ds_path = path
        
        self.gene_index = pl.read_csv(gene_index_path, sep = '\t')
        
        dep_map_data = pl.read_csv(path, sep = '\t')
        
        dep_map_data = dep_map_data.melt(id_vars=['gene_id'],
                                         variable_name='cell_line',
                                         value_name='CERES_score',
                                         )
        
        dep_map_data = dep_map_data.with_column(pl.col('cell_line').cast(pl.Categorical))

        self.cell_line_index = dep_map_data.select(pl.col('cell_line').unique())
        self.cell_line_index = self.cell_line_index.with_row_count(name = 'cell_line_idx')
        
        dep_map_data = dep_map_data.join(self.gene_index, on='gene_id')
        dep_map_data = dep_map_data.join(self.cell_line_index, on='cell_line')
        
        # normalize by deviding by global std
        dep_map_data = dep_map_data.with_column(pl.col('CERES_score') / pl.col('CERES_score').std())
        
        # shuffle and split in train / test 
        dep_map_data = dep_map_data.sample(frac = 1, shuffle = True)
        
        self.train_table = dep_map_data[:int(len(dep_map_data) * 0.9)]
        self.test_table = dep_map_data[int(len(dep_map_data) * 0.9):]
        
                
    def compute_sample_init_pca(self, dim = 32, return_loadings = False, scale_pcs = True):
        
        dep_map_data = pl.read_csv(self.ds_path, sep = '\t')
        pca = PCA(dim)
        
        ds = dep_map_data.select(self.cell_line_index.get_column('cell_line').to_list()).to_numpy().transpose()
        ds = ds / ds.std()
        
        pcs = pca.fit_transform(ds)
        
        if scale_pcs is True:
            pcs = pcs / pcs.std(axis = 0)
        
        if return_loadings is True:
            return pcs, pca.components_.transpose(), dep_map_data.join(self.gene_index, on = 'gene_id', how='left').gene_idx.to_numpy()
        else: 
            return pcs
    
    
    def genes_in_assay(self):
        return self.train_table.get_column("gene_id").unique()
    
    
    def get_num_batches_per_epoch(self, batchsize):
        return int(np.ceil(len(self.train_table) / batchsize))
    
        
    def get_batches(self, batchsize, subset = 'train'):
        
        if subset == "train":
            table = self.train_table
        if subset == "test":
            table = self.test_table
        
        table = table.sample(frac = 1, shuffle = True)
        # compute number of batches
        n_batches = int(np.ceil(len(table) / batchsize))
        
        start_idx = 0
        
        for batch_idx in range(n_batches):
            if batch_idx < (n_batches - 1):
                s = slice(start_idx, start_idx + batchsize)
                start_idx += batchsize
            else:
                s = slice(start_idx, len(table))

            subset_table = table[s]
            g_idx = torch.from_numpy(subset_table.get_column("gene_idx").to_numpy().astype('int64'))
            s_idx = torch.from_numpy(subset_table.get_column("cell_line_idx").to_numpy().astype('int64'))
            score = torch.from_numpy(subset_table.get_column("CERES_score").to_numpy().astype('float32'))
            yield g_idx, s_idx, score
    
    

class ProtT5Dataset(Dataset):
    def __init__(self,
                 path = project_path + '/embedding/individual/prot_t5_128_embedding.tsv',
                 gene_index_path = project_path + '/processed_data/omics_index/gene_index.tsv'):
        
        self.ds_path = path
        
        self.gene_index = pl.read_csv(gene_index_path, sep = '\t')
        
        prot_t5 = pl.read_csv(path, sep = '\t')
        
        prot_t5 = prot_t5.join(self.gene_index, on='gene_id')
        self.dataset = prot_t5
     
        
    def __len__(self): 
        return len(self.dataset)
    
    def __getitem__(self, idx):
        
        subset_table = self.dataset[idx]
        g_idx = torch.from_numpy(subset_table.gene_idx.to_numpy().astype('int64'))
        emb = torch.from_numpy(subset_table.drop(['gene_id', 'gene_idx']).to_numpy().astype('float32'))
        assert emb.shape[1] == 128
        
        return g_idx, emb
    
    
 
    
class FasterGTExDataset(Dataset):
    def __init__(self, 
                 path = project_path + '/processed_data/gtex/log10_tpm_subset.tsv',
                 gene_index_path = project_path + '/processed_data/omics_index/gene_index.tsv',
                 ):
        
        self.ds_path = path
        
        gtex_data = pl.read_csv(path, sep = '\t')
        
        gtex_data = gtex_data.melt(id_vars=['gene_id'],
                           variable_name='sample_id',
                           value_name='TPM_norm_score',
                           )
        
        gtex_data = gtex_data.with_column(pl.col('sample_id').cast(pl.Categorical))
        
        self.gtex_sample_index = gtex_data.select(
                pl.col("sample_id").unique()
        )
        self.gtex_sample_index = self.gtex_sample_index.with_row_count(name = "sample_idx")
        
        self.gene_index = pl.read_csv(gene_index_path, sep = '\t')
        
        gtex_data = gtex_data.join(self.gtex_sample_index, on = 'sample_id')
        gtex_data = gtex_data.join(self.gene_index, on = 'gene_id')
        
        
        # center gtex
        gtex_data = gtex_data.with_column(
            pl.col('TPM_norm_score') - pl.col('TPM_norm_score').mean().over('gene_id')
        )
        gtex_data = gtex_data.with_column(
            pl.col('TPM_norm_score') / pl.col('TPM_norm_score').std()
        )
        
        self.dataset = gtex_data
        
        # shuffle and split in train / test 
        self.dataset = self.dataset.sample(frac = 1, shuffle = True)
        
        self.train_table = self.dataset[:int(len(self.dataset) * 0.9)]
        self.test_table = self.dataset[int(len(self.dataset) * 0.9):]
        
        
    def compute_sample_init_pca(self, dim = 64, return_loadings=False):
        
        gtex_data = pl.read_csv(self.ds_path, sep = '\t')
        pca = PCA(dim)
        
        # center genes and divide by global std.
        ds = gtex_data.select(self.gtex_sample_index.get_column('sample_id').to_list()).to_numpy()
        ds = ds - ds.mean(axis = 1).reshape(-1, 1)
        ds = ds / ds.std()
        
        pcs = pca.fit_transform(ds.transpose())
        
        if return_loadings is True:
            return pcs, pca.components_.transpose(), gtex_data.join(self.gene_index, on = 'gene_id', how='left').gene_idx.to_numpy()
        else: 
            return pcs
    
    def get_num_batches_per_epoch(self, batchsize):
        return int(np.ceil(len(self.train_table) / batchsize))
    
    
    def genes_in_assay(self):
        return self.train_table.get_column("gene_id").unique()
    
    
    def get_batches(self, batchsize, subset = 'train'):
        
        if subset == "train":
            table = self.train_table
        if subset == "test":
            table = self.test_table
        
        table = table.sample(frac = 1, shuffle = True)
        # compute number of batches
        n_batches = int(np.ceil(len(table) / batchsize))
        
        start_idx = 0
        
        for batch_idx in range(n_batches):
            if batch_idx < (n_batches - 1):
                s = slice(start_idx, start_idx + batchsize)
                start_idx += batchsize
            else:
                s = slice(start_idx, len(table))

            subset_table = table[s]
            g_idx = torch.from_numpy(subset_table.get_column("gene_idx").to_numpy().astype('int64'))
            s_idx = torch.from_numpy(subset_table.get_column("sample_idx").to_numpy().astype('int64'))
            score = torch.from_numpy(subset_table.get_column("TPM_norm_score").to_numpy().astype('float32'))
            yield g_idx, s_idx, score
            
            
            
class ProtT5Dataset(Dataset):
    def __init__(self, 
                 path = project_path + '/embedding/individual/prot_t5_embedding.tsv',
                 gene_index_path = project_path + '/processed_data/omics_index/gene_index.tsv',
                 ):
        
        self.ds_path = path
        
        data = pl.read_csv(path, sep = '\t')
        
        data = data.melt(id_vars=['gene_id'],
                         variable_name='sample_id',
                         value_name='score',
                          )
        
        data = data.with_column(pl.col('sample_id').cast(pl.Categorical))
        
        self.sample_index = data.select(
                pl.col("sample_id").unique()
        )
        self.sample_index = self.sample_index.with_row_count(name = "sample_idx")
        
        self.gene_index = pl.read_csv(gene_index_path, sep = '\t')
        
        data = data.join(self.sample_index, on = 'sample_id')
        data = data.join(self.gene_index, on = 'gene_id')
        
        self.dataset = data
        
        # shuffle and split in train / test 
        self.dataset = self.dataset.sample(frac = 1, shuffle = True)
        
        self.train_table = self.dataset[:int(len(self.dataset) * 0.9)]
        self.test_table = self.dataset[int(len(self.dataset) * 0.9):]
        
        
    def compute_sample_init_pca(self, dim = 128, return_loadings=False):
        
        data = pl.read_csv(self.ds_path, sep = '\t')
        pca = PCA(dim)
        
        ds = data.select(self.sample_index.get_column('sample_id').to_list()).to_numpy()
        
        pcs = pca.fit_transform(ds.transpose())
        
        if return_loadings is True:
            return pcs, pca.components_.transpose(), data.join(self.gene_index, on = 'gene_id', how='left').gene_idx.to_numpy()
        else: 
            return pcs
    
    def get_num_batches_per_epoch(self, batchsize):
        return int(np.ceil(len(self.train_table) / batchsize))
    
    
    def genes_in_assay(self):
        return self.train_table.get_column("gene_id").unique()
    
    
    def get_batches(self, batchsize, subset = 'train'):
        
        if subset == "train":
            table = self.train_table
        if subset == "test":
            table = self.test_table
        
        table = table.sample(frac = 1, shuffle = True)
        # compute number of batches
        n_batches = int(np.ceil(len(table) / batchsize))
        
        start_idx = 0
        
        for batch_idx in range(n_batches):
            if batch_idx < (n_batches - 1):
                s = slice(start_idx, start_idx + batchsize)
                start_idx += batchsize
            else:
                s = slice(start_idx, len(table))

            subset_table = table[s]
            g_idx = torch.from_numpy(subset_table.get_column("gene_idx").to_numpy().astype('int64'))
            s_idx = torch.from_numpy(subset_table.get_column("sample_idx").to_numpy().astype('int64'))
            score = torch.from_numpy(subset_table.get_column("score").to_numpy().astype('float32'))
            yield g_idx, s_idx, score
            
  