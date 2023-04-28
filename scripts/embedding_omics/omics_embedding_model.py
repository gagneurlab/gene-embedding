#!/usr/bin/env python
# coding: utf-8
# %%

# %%

import optuna

import gc
import os
import torch
import numpy as np
import polars as pl
import umap 

import scipy.stats

import random

from torch.utils.data import Dataset, DataLoader, BatchSampler, RandomSampler, random_split
import plotnine as p9

import torch_geometric
from torch_geometric.utils import to_undirected
from scvi.distributions import NegativeBinomial, Poisson, ZeroInflatedNegativeBinomial

from sklearn.metrics import r2_score, roc_auc_score, average_precision_score
import wandb
import sys


sys.path.append('../../pytorch-revgrad/src/')
sys.path.append(os.getcwd())

from pytorch_revgrad import RevGrad


device = 'cuda'

NUM_WORKERS = 0


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--models_to_train", type=str, nargs="*", help="models which should be trained")
parser.add_argument("--EPOCHS", type=int, default = 500)
parser.add_argument("--LEARNING_RATE", type=float, default = 1e-3)
parser.add_argument("--GENE_EMB_SIZE", type=int, default = 265)
parser.add_argument("--DEP_MAP_EMB_SIZE", type=int, default = 64)
parser.add_argument("--DEP_MAP_LAMBDA", type=float, default = 1)
parser.add_argument("--GTEX_EMB_SIZE", type=int, default = 256)
parser.add_argument("--PROT_T5_EMB_SIZE", type=int, default = 128)
parser.add_argument("--NUM_LAYERS", type=int, default = 3)
parser.add_argument("--domain_adapt", type=bool, default = True)
parser.add_argument("--DOMAIN_ADAPT_LAMBDA", type=float, default = 1.0)
parser.add_argument("--DOMAIN_ADAPT_LAYERS", type=int, default = 2)
parser.add_argument("--DOMAIN_ADAPT_LOSS", type=str, default = 'dot_prod')
parser.add_argument("--DOMAIN_ADAPT_RSAMPLE", type=str, default = "no")
parser.add_argument("--ASSAY_VAR_INIT_FACTOR", type=float, default = 2.)
parser.add_argument("--WEIGHT_TS_LOSS", type=bool, default = False)
parser.add_argument("--OPT", type=str, default = 'adam')
parser.add_argument("--OUTPATH", type=str, default = None)


args = parser.parse_args()

# %%
models_to_train  =  args.models_to_train
EPOCHS = args.EPOCHS
LEARNING_RATE  =  args.LEARNING_RATE
GENE_EMB_SIZE  =  args.GENE_EMB_SIZE
DEP_MAP_EMB_SIZE  =  args.DEP_MAP_EMB_SIZE
DEP_MAP_LAMBDA = args.DEP_MAP_LAMBDA
GTEX_EMB_SIZE  =  args.GTEX_EMB_SIZE
PROT_T5_EMB_SIZE = args.PROT_T5_EMB_SIZE
NUM_LAYERS = args.NUM_LAYERS
domain_adapt  =  args.domain_adapt
DOMAIN_ADAPT_LAMBDA = args.DOMAIN_ADAPT_LAMBDA
DOMAIN_ADAPT_LAYERS = args.DOMAIN_ADAPT_LAYERS
DOMAIN_ADAPT_RSAMPLE = args.DOMAIN_ADAPT_RSAMPLE
DOMAIN_ADAPT_LOSS = args.DOMAIN_ADAPT_LOSS
ASSAY_VAR_INIT_FACTOR = args.ASSAY_VAR_INIT_FACTOR
WEIGHT_TS_LOSS = args.WEIGHT_TS_LOSS
OPT = args.OPT
OUTPATH = args.OUTPATH


print("Models to train:")
print(models_to_train)

models_to_train_short = models_to_train.copy()
    
print(models_to_train)

from omics_datasets import *
from  torch.distributions.beta import Beta
from torch.distributions.normal import Normal
from torch.distributions.kl import kl_divergence
import torch.nn as nn


class GeneEmbedding(torch.nn.Module):
    def __init__(self, n_genes, emb_dim, gene_list, dir_concentration = 0.1, emb_init = None):
        super(GeneEmbedding, self).__init__()
        
        self.n_genes = n_genes
        self.emb_dim = emb_dim
        self.genes = gene_list
        
    
        if emb_init is not None:
            self.emb_mu = nn.Embedding.from_pretrained(emb_init, freeze = False)
        else:
            # init as [0, 1], as the dirichlet prior is in [0, 1]
            self.emb_mu = nn.Embedding.from_pretrained(
                torch.rand((n_genes, emb_dim)), freeze = False
            )
    
        self.emb_log_sigma = nn.Embedding.from_pretrained(
            torch.full((n_genes, emb_dim), np.log(0.5)), freeze = False
        )
    
    
    def get_emb_table(self):
        emb_df = pd.DataFrame(self.emb_mu.weight.detach().cpu(), 
                              index = self.genes,
                              columns = [f'FACT_EMB_{i}' for i in range(self.emb_dim)]
                              )
        emb_df.index.name = 'gene_id'
        return emb_df
    
    
    def get_log_sigma_table(self):
        emb_df = pd.DataFrame(self.emb_log_sigma.weight.detach().cpu(), 
                              index = self.genes,
                              columns = [f'FACT_EMB_{i}' for i in range(self.emb_dim)]
                              )
        emb_df.index.name = 'gene_id'
        return emb_df
    
    
    def get_shape(self):
        return (self.n_genes, self.emb_dim)
    

    
class AssayModel(torch.nn.Module):
    def __init__(self, gene_emb, sigma_init = np.log(0.5)):
        super().__init__()
        
        self.gene_emb = gene_emb
        
        self.assay_log_sigma = nn.Embedding.from_pretrained(
            torch.full(gene_emb.get_shape(), sigma_init), freeze = False
        )
        
        self.pass_assay_log_sigma = False
        
        
    def reparameterization(self, mean, sd):
        epsilon = torch.randn_like(sd)    # sampling epsilon        
        z = mean + sd * epsilon           # reparameterization trick
        return z


    # get mean log_sigma of the current assay without gradient. 
    def get_assay_log_sigma_mean(self, idx):
        with torch.no_grad():
            assay_log_sigma_mean = self.assay_log_sigma.weight.mean(axis = 1, keepdim = True)
            assay_log_sigma_mean = assay_log_sigma_mean[idx,:]
        return assay_log_sigma_mean
    
    def get_emb(self, idx):
        emb_mu = self.gene_emb.emb_mu.weight
        emb_sigma = self.gene_emb.emb_log_sigma.weight.exp()
        assay_sigma = self.assay_log_sigma.weight.exp()
        
                # replace the sampling twice by sampling only once using mu = mu_1 + mu_2 and sig^2 = sig_1^2 + sig_2^2
        sigma = torch.sqrt(emb_sigma.pow(2) + assay_sigma.pow(2))
        z_emb_assay = self.reparameterization(emb_mu[idx, :], sigma[idx, :])
        
        # add assay log sigma to embedding without gradient.
        if self.pass_assay_log_sigma is True:
            assay_log_sigma_mean = self.get_assay_log_sigma_mean(idx)
            z_emb_assay = torch.cat((z_emb_assay, assay_log_sigma_mean), dim = 1)
        return z_emb_assay
    
    
    def get_assay_log_sigma_table(self):
        emb_df = pd.DataFrame(self.assay_log_sigma.weight.detach().cpu(), 
                              index = self.gene_emb.genes,
                              columns = [f'FACT_EMB_{i}' for i in range(self.gene_emb.emb_dim)]
                              )
        emb_df.index.name = 'gene_id'
        return emb_df
    
### -------------------------------- Tabular Data Model and Loss ------------------------------------------ ### 


class TabluarDataModel(AssayModel):
    def __init__(self,  gene_emb, n_samples, emb_dim, sigma_init, pass_assay_log_sigma = False, emb_init = None):
        super().__init__(gene_emb, sigma_init)
            
        self.pass_assay_log_sigma = pass_assay_log_sigma
        
        if emb_init is not None:
            self.sample_emb = nn.Embedding.from_pretrained(emb_init, freeze=False)
            #self.register_parameter("sample_emb",nn.Embedding.from_pretrained(emb_init, freeze=False))
        else:
            self.sample_emb = nn.Embedding(n_samples, emb_dim)
        
        joint_emb_dim = emb_dim + self.gene_emb.emb_dim
        var_model_dim = self.gene_emb.emb_dim
        
        if pass_assay_log_sigma is True:
            joint_emb_dim += 1
            var_model_dim += 1

        self.model = nn.Sequential()
        
        for i in range(NUM_LAYERS - 1):
            self.model.add_module(f"layer_{i}", nn.Linear(joint_emb_dim, joint_emb_dim))
            self.model.add_module(f"relu_{i}", nn.ReLU())
        
        self.model.add_module(f"layer_{NUM_LAYERS - 1}", nn.Linear(joint_emb_dim, 1))
    
        
        self.var_model = nn.Sequential()
        
        for i in range(NUM_LAYERS - 1):
            self.var_model.add_module(f"layer_{i}", nn.Linear(var_model_dim, var_model_dim))
            self.var_model.add_module(f"relu_{i}", nn.ReLU())
        
        self.var_model.add_module(f"layer_{NUM_LAYERS - 1}", nn.Linear(var_model_dim, 1))
        # Init bias of last layer with small variance.
        nn.init.constant_(self.var_model[-1].bias, -2)
    
    #@torch.jit.script
    def forward(self, gene_index, sample_index):
        
        gene_emb_batch = self.get_emb(gene_index)
        sample_emb_batch = self.sample_emb(sample_index)
        
        joint_emb = torch.cat((gene_emb_batch, sample_emb_batch), dim = 1)
        
        pred = self.model(joint_emb)
        log_var = self.var_model(gene_emb_batch)
        # clip to avoid exp of log var to become inf or 0.
        log_var = log_var.clip(min=-40, max=40)
        
        return pred.squeeze(), log_var.squeeze()

    

# Domain adaptation model
class PredictAssay(torch.nn.Module):
    def __init__(self, gene_emb, num_assays):
        super().__init__()
        
        self.gene_emb = gene_emb
        in_dim = self.gene_emb.emb_mu.weight.shape[1]
        
        self.model = nn.Sequential(RevGrad())
        
        for i in range(DOMAIN_ADAPT_LAYERS - 1):
            self.model.add_module(f"layer_{i}", nn.Linear(in_dim, in_dim))
            self.model.add_module(f"reul_{i}", nn.ReLU())
        
        self.model.add_module(f"layer_{DOMAIN_ADAPT_LAYERS}", nn.Linear(in_dim, num_assays))
        
        self.random_sample = DOMAIN_ADAPT_RSAMPLE
        
        
    def reparameterization(self, mean, sd):
        epsilon = torch.randn_like(sd)    # sampling epsilon        
        z = mean + sd * epsilon           # reparameterization trick
        return z
    
    
    def forward(self, assay_sig):
        
        if self.random_sample == "gene_sig":
            mu = self.gene_emb.emb_mu.weight
            sigma = self.gene_emb.emb_log_sigma.weight.exp()
            
            emb = self.reparameterization(mu, sigma)
            
        elif self.random_sample == "gene_assay_sig":
            mu = self.gene_emb.emb_mu.weight
            sigma = self.gene_emb.emb_log_sigma.weight.exp() + assay_sig
            
            emb = self.reparameterization(mu, sigma)
 
        else:
            emb = self.gene_emb.emb_mu.weight
        
        
        return self.model(emb)

# %%

class LimitedBCE(nn.Module):
    def __init__(self):
        super().__init__()
        
        self.loss = nn.BCEWithLogitsLoss()
        
        
    def forward(self, logit_pred, lables):
        naive_pred = lables.sum() / lables.numel()
        naive_pred = torch.logit(naive_pred)
        naive_loss = self.loss(naive_pred.expand_as(lables), lables)
        loss = self.loss(logit_pred, lables)
        print(naive_loss)
        
        return min(naive_loss, loss)
    
    
class NegAbsDotProd(nn.Module):
    def __init__(self):
        super().__init__()
        
    def forward(self, logit_pred, lables):
        lables = lables * 2 - 1
        
        #print((logit_pred * lables).abs())
        abs_dot_prod = (logit_pred * lables).abs()
        
        return - abs_dot_prod.mean()
    

class NegAbsDotProd2(nn.Module):
    def __init__(self):
        super().__init__()
        
    def forward(self, logit_pred, lables):
        # cast from [0, 1] to [-1, 1]
        lables = lables * 2 - 1
        
        #print((logit_pred * lables).abs())
        dot_prod = (logit_pred * lables)
        abs_dot_prod = dot_prod.mean(axis=1).abs()
        
        return - abs_dot_prod.mean()


def compute_string_average_precision(emb):
    dist_frame = pl.read_csv('...path_to_project...' + '/processed_data/string_hyperopt/string.tsv',
                    sep = '\t')
    
    dl = DataLoader(np.arange(len(dist_frame)), batch_size = 100_000)
    dist = []
    for b in dl:
        temp = dist_frame[np.array(b)]
        dist.append(np.sqrt(((emb.reindex(temp.get_column('gene_id1')).reset_index(drop=True) - 
                              emb.reindex(temp.get_column('gene_id2')).reset_index(drop=True))**2).sum(axis = 1)))
        
    dist_frame = dist_frame.with_column(pl.Series(np.concatenate(dist)).alias('pdist'))
    return average_precision_score(dist_frame.get_column('score'), -dist_frame.get_column('pdist'))


device = 'cuda'


# %%


NUM_WORKERS = 0


# %%


datasets = {}
for m in models_to_train:
    
    if m == "depMap":
        datasets[m] = FasterDepMapDataset()
        
    elif m == "gtex":
        datasets[m] = FasterGTExDataset()
    
    elif m == "protT5":
        datasets[m] = ProtT5Dataset()
        
    else:
        print(f"Dataset {m} not implemented")
        break


# %%
if domain_adapt is True:
    target = datasets[next(iter(datasets))].gene_index.select('gene_id')
    
    for m in models_to_train:
        target = target.with_column(pl.col('gene_id')
                                    .is_in(
                    datasets[m].genes_in_assay()
                                          )
                                    .alias(f"in_{m}")
                                     )
 
    domain_adapt_target = torch.from_numpy(target.drop('gene_id').to_numpy().astype('float32'))
        
        
    

# %%
domain_adapt_target_device = domain_adapt_target.to(device)

# %%


datasets


print("Datasets loaded!!")


batchsizes = {"depMap": int(1e6),
              "protT5": int(1e6),
              "gtex": int(1e6),
             }



###  ----------------------- Variational Loss functions ----------------------------------------


class BaseLoss(nn.Module):
    def __init__(self, num_assays, assay_prior_log_sigma, mu, rho):
        super(BaseLoss, self).__init__()
        
        self.num_assays = torch.as_tensor(num_assays, device=device)
        self.assay_ls = torch.as_tensor(assay_prior_log_sigma, device=device)
        self.emb_log_sigma_prior = torch.as_tensor(0., device=device)
        
        
    def compute_kl(self, model):
        
        loss_KLD_emb = kl_divergence(Normal(model.gene_emb.emb_mu.weight, model.gene_emb.emb_log_sigma.weight.exp()),
                                     Normal(0, self.emb_log_sigma_prior.exp())
                                     ).mean()
        loss_KLD_assay = kl_divergence(Normal(0, model.assay_log_sigma.weight.exp()),
                                       Normal(0, self.assay_ls.exp())
                                       ).mean()
        return loss_KLD_assay + loss_KLD_emb / self.num_assays 



class MatLoss(BaseLoss):
    def __init__(self, num_assays, assay_prior_log_sigma, mu, rho):
        super(MatLoss, self).__init__(num_assays, assay_prior_log_sigma, mu, rho)
        
        self.mse_loss_mean = nn.MSELoss(reduction="mean")
        self.gnll = nn.GaussianNLLLoss(full=True, reduction='mean')
    
        
    def forward(self, pred_mu, pred_log_var, target, model):
           
        loss_NLL = self.gnll(pred_mu, target, torch.exp(pred_log_var))
        kl_term = self.compute_kl(model)
        
        return loss_NLL + kl_term, self.mse_loss_mean(pred_mu, target), loss_NLL #, loss_KLD_emb + loss_KLD_assay




### ------------------------ train functions ---------------------------------------

# %%
def train_epoch(model, opt, loss_func, d_adapt_model, d_adapt_loss, dataset, dataset_name, lambda_model):
    def train(epoch):
        epoch_elbo = []
        epoch_mse = []
        epoch_nll = []
        epoch_val_elbo = []
        epoch_val_mse = []
        epoch_val_nll = []
        epoch_domain_adapt = []
        
        first_batch = True
        #model.to(device)
        for g_idx, s_idx, target in dataset.get_batches(int(1e6), 'train'):
            print(g_idx)
            print(s_idx)
            print(target)
            
            print(g_idx.shape)
            print(s_idx.shape)
            print(target.shape)
            
            
            opt.zero_grad()
            mu, log_var = model(g_idx.to(device, non_blocking=True), s_idx.to(device, non_blocking=True))
            target = target.to(device, non_blocking=True)
            elbo, mse, nll = loss_func(mu, 
                                       log_var, 
                                       target,
                                       model
                                       )

            
            epoch_elbo.append(elbo.item())
            epoch_mse.append(mse.item())
            epoch_nll.append(nll.item())
            
            # compute domain adaptation loss
            domain_adapt_pred = d_adapt_model(model.assay_log_sigma.weight.exp())
            domain_adapt_loss = d_adapt_loss(domain_adapt_pred, domain_adapt_target_device)
            epoch_domain_adapt.append(domain_adapt_loss.item())
            
            # sum losses
            loss = lambda_model * elbo + domain_adapt_loss * DOMAIN_ADAPT_LAMBDA
            loss.backward()
            opt.step()
            
            if first_batch is True:
                wandb.log({f'{dataset_name}_grad_mean_emb': wandb.Histogram(model.gene_emb.emb_mu.weight.grad.detach().cpu())}, step = epoch, commit = False)
                first_batch = False
                        
        opt.zero_grad(set_to_none=True)
        torch.cuda.empty_cache()
        
        val_pred = []
        val_target = []
        with torch.no_grad():
             for g_idx, s_idx, target in dataset.get_batches(int(1e6), 'test'):
                    mu, log_var = model(g_idx.to(device, non_blocking=True), s_idx.to(device, non_blocking=True))
                    target = target.to(device, non_blocking=True)
                    elbo, mse, nll = loss_func(mu, 
                                               log_var, 
                                               target,
                                               model
                                               )

                    val_pred.append(mu.detach().cpu().numpy())
                    val_target.append(target.detach().cpu().numpy())
                    epoch_val_elbo.append(elbo.detach().item())
                    epoch_val_mse.append(mse.detach().item())
                    epoch_val_nll.append(nll.detach().item())

        ## compute r2
        val_pred = np.concatenate(val_pred)
        val_target = np.concatenate(val_target)

        r2 = r2_score(val_target, val_pred)
        
        print(f"Epoch {epoch} done, {dataset_name} Train ELBO: {np.round(np.mean(epoch_elbo),3)}, Val ELBO {np.round(np.mean(epoch_val_elbo),3)}, Train MSE {np.round(np.mean(epoch_mse),3)}, Val MSE {np.round(np.mean(epoch_val_mse),3)}")
        
        wandb.log({f"{dataset_name}_ELBO": np.mean(epoch_elbo), 
                   f"{dataset_name}_val_ELBO": np.mean(epoch_val_elbo),
                   f"{dataset_name}_NLL": np.mean(epoch_nll),
                   f"{dataset_name}_val_NLL": np.mean(epoch_val_nll),
                   f"{dataset_name}_MSE": np.mean(epoch_mse),
                   f"{dataset_name}_val_MSE": np.mean(epoch_val_mse),
                   f"{dataset_name}_domain_adapt": np.mean(epoch_domain_adapt),
                  }, 
                  step = epoch,
                  commit = False
                  )
        
        #model = model.cpu()
        return r2
    return train
        

    
##  ------------------- init models ---------------------


def objective(trial):
    
    if SAME_ASSAY_PRIOR is True:

        l_prior = 0

        assay_prior_log_sigma = {
            "mu": 0,
            "rho": 0,
            #"gene_emb_var_prior": gene_emb_var_prior,
            "depMap": l_prior,
            "gtex": l_prior,
            "protT5": l_prior,
        }
    else:
        l_depMap = trial.suggest_float('l_depMap', -3, 3)
        l_gtex = trial.suggest_float('l_gtex', -3, 3)
        l_protT5 = trial.suggest_float('l_protT5', -3, 3)

        assay_prior_log_sigma = {
            "mu": 0,
            "rho": 0,
            #"gene_emb_var_prior": gene_emb_var_prior,
            "depMap": l_depMap,
            "gtex": l_gtex,
            "protT5": l_protT5,
        }




    n_genes = len(next(iter(datasets.values())).gene_index)
    emb_dim = GENE_EMB_SIZE
    gene_emb = GeneEmbedding(n_genes, 
                             emb_dim,  
                             next(iter(datasets.values())).gene_index.get_column('gene_id').to_list(), 
                             emb_init = torch.rand(n_genes, emb_dim) - 0.5,
                             #dir_concentration = assay_prior_log_sigma['d_concentration']
                            )


    models = {}
    for m in models_to_train:


        if m  == "depMap":
            models[m] = TabluarDataModel(gene_emb, 
                                    n_samples = len(datasets['depMap'].cell_line_index), 
                                    emb_dim = DEP_MAP_EMB_SIZE,
                                    # emb_init = torch.from_numpy(depMap_sample_init.astype('float32'))
                                    sigma_init = assay_prior_log_sigma['depMap'] - np.log(ASSAY_VAR_INIT_FACTOR),
                                    emb_init = torch.from_numpy(datasets['depMap'].compute_sample_init_pca(DEP_MAP_EMB_SIZE).astype('float32')),
                                    pass_assay_log_sigma = PASS_ASSAY_LOG_SIGMA,
                                    )

        if m  == "gtex":
            models[m] = TabluarDataModel(gene_emb, 
                                  n_samples = len(datasets['gtex'].gtex_sample_index), 
                                  emb_dim = GTEX_EMB_SIZE,
                                  # emb_init = torch.from_numpy(gtex_sample_init.astype('float32'))
                                    sigma_init = assay_prior_log_sigma['gtex'] - np.log(ASSAY_VAR_INIT_FACTOR),
                                    emb_init = torch.from_numpy(datasets['gtex'].compute_sample_init_pca(GTEX_EMB_SIZE).astype('float32')),
                                    pass_assay_log_sigma = PASS_ASSAY_LOG_SIGMA,
                                  )

        if m == "protT5":
            models[m] = TabluarDataModel(gene_emb, 
                                    n_samples = len(datasets['protT5'].sample_index), 
                                    emb_dim = PROT_T5_EMB_SIZE,
                                  # emb_init = torch.from_numpy(gtex_sample_init.astype('float32'))
                                    sigma_init = assay_prior_log_sigma['depMap'] - np.log(ASSAY_VAR_INIT_FACTOR),
                                    emb_init = torch.from_numpy(datasets['protT5'].compute_sample_init_pca(PROT_T5_EMB_SIZE).astype('float32')),
                                    pass_assay_log_sigma = PASS_ASSAY_LOG_SIGMA,
                                    )

            
      
    if domain_adapt is True:
        models['domainAdapt'] = PredictAssay(gene_emb, len(models_to_train))


    losses = {"depMap": MatLoss(len(models_to_train_short), 
                                assay_prior_log_sigma['depMap'], 
                                assay_prior_log_sigma['mu'], 
                                assay_prior_log_sigma['rho']),
              "gtex": MatLoss(len(models_to_train_short), 
                              assay_prior_log_sigma['gtex'], 
                              assay_prior_log_sigma['mu'], 
                              assay_prior_log_sigma['rho']),
              "protT5": MatLoss(len(models_to_train_short), 
                                assay_prior_log_sigma['protT5'], 
                                assay_prior_log_sigma['mu'], 
                                assay_prior_log_sigma['rho']),
             }

    print(DOMAIN_ADAPT_LOSS)
    if DOMAIN_ADAPT_LOSS == "bce":
        losses["domainAdapt"] = nn.BCEWithLogitsLoss()
    elif DOMAIN_ADAPT_LOSS == "dot_prod":
        losses["domainAdapt"] = NegAbsDotProd() 
    else:
        print("Domain adaptation loss not implemented")


    print("Model init done!!")



    num_batches = {}
    for m in models_to_train:
        num_batches[m] = datasets[m].get_num_batches_per_epoch(batchsizes[m])
        print( datasets[m].get_num_batches_per_epoch(batchsizes[m]))



    wandb.init(project=f"OMICS_EMBEDDING",
               entity=".....set your entity ....", 
               reinit = True, 
               config = {
                   "models_to_train":  models_to_train_short,
                   "DEP_MAP_EMB_SIZE": DEP_MAP_EMB_SIZE,
                   "GTEX_EMB_SIZE": GTEX_EMB_SIZE,
                   "PROT_T5_EMB_SIZE": PROT_T5_EMB_SIZE,
                   "GENE_EMB_SIZE": GENE_EMB_SIZE,
                   "domain_adapt":domain_adapt,
                   "NUM_LAYERS": NUM_LAYERS,
                   "LEARNING_RATE": LEARNING_RATE,
                   "DOMAIN_ADAPT_LAMBDA": DOMAIN_ADAPT_LAMBDA,
                   "DOMAIN_ADAPT_LAYERS": DOMAIN_ADAPT_LAYERS,
                   "DOMAIN_ADAPT_LOSS": DOMAIN_ADAPT_LOSS,
                   "ASSAY_VAR_INIT_FACTOR": ASSAY_VAR_INIT_FACTOR,
                   "EPOCHS": EPOCHS,
                   "OPT": OPT,
                   "assay_prior_log_sigma": assay_prior_log_sigma,
                   "SLURM_JOB_ID": os.environ['SLURM_JOB_ID'],
                   "num_batches": num_batches,
               })

    for m in models.values():
        m.to(device)

    # %%
    opt = {}
    for m in models:



        if OPT == "ranger":
            if m == "domainAdapt":
                num_batches_per_epoch = 1
            else:
                num_batches_per_epoch = datasets[m].get_num_batches_per_epoch(batchsizes[m])
            opt[m] = Ranger21(models[m].parameters(), 
                              lr = LEARNING_RATE, 
                              num_epochs = EPOCHS, 
                              num_batches_per_epoch=num_batches_per_epoch, 
                              )
        if OPT == "adam":
            opt[m] = torch.optim.Adam(models[m].parameters(), 
                                      lr = LEARNING_RATE,
                                      )

        if OPT == "adam_lr_adj":

            # create lists of params.
            emb_param_names = ['gene_emb.emb_mu.weight', 'gene_emb.emb_log_sigma.weight', 'gene_emb.concentration']
            emb_params = [kv[1] for kv in models[m].named_parameters() if kv[0] in emb_param_names]
            model_params = [kv[1] for kv in models[m].named_parameters() if kv[0] not in emb_param_names]

            opt[m] = torch.optim.Adam([
                                            {'params': emb_params, "lr": LEARNING_RATE / len(models_to_train)},
                                            {'params': model_params}
                                      ],
                                      lr = LEARNING_RATE,
                                      )

        if OPT == "adam_lr_bs_sqrt":

            # create lists of params.
            emb_param_names = ['gene_emb.emb_mu.weight', 'gene_emb.emb_log_sigma.weight']
            emb_params = [kv[1] for kv in models[m].named_parameters() if kv[0] in emb_param_names]
            model_params = [kv[1] for kv in models[m].named_parameters() if kv[0] not in emb_param_names]

            # get number of batches for each dataset
            if m.startswith("tabulaSapiensTissue") & (WEIGHT_TS_LOSS is False):    
                bs = num_batches["tabulaSapiensTissue"]
            elif m == "domainAdapt":
                bs = 1
            else:
                bs = num_batches[m]

            opt[m] = torch.optim.Adam([
                                            {'params': emb_params, "lr": LEARNING_RATE / np.sqrt(bs)},
                                            {'params': model_params}
                                      ],
                                      lr = LEARNING_RATE,
                                      )


    print(opt)

    # %%
    if OUTPATH is not None:
        exp_path = OUTPATH
    else:
        exp_path = f"...path_to_embeddings..../{wandb.run.name}_{os.environ['SLURM_JOB_ID']}/"

    # %%
    from pathlib import Path

    # %%
    Path(exp_path).mkdir(parents=True, exist_ok=True)


    # %%
    functions = []

    for m in models_to_train:
        if m == "depMap":
            functions.append(
                {"name": m,
                 "train_func": train_epoch(models['depMap'], opt['depMap'], losses['depMap'], models["domainAdapt"], losses["domainAdapt"], datasets['depMap'], "depMap", 1.)
                }
            )
        elif m == "gtex":
            functions.append(
                {"name": m,
                 "train_func": train_epoch(models['gtex'], opt['gtex'], losses['gtex'], models["domainAdapt"], losses["domainAdapt"], datasets['gtex'], "gtex", 1.0)
                }
            ) 
        elif m == "protT5":
            functions.append(
                {"name": m,
                 "train_func": train_epoch(models['protT5'], opt['protT5'], losses['protT5'], models["domainAdapt"], losses["domainAdapt"], datasets['protT5'], "protT5", 1.0)
                }
            )
        else:
            print(f"Train function for {m} not implemented")
            break

    print(functions)


    torch.cuda.empty_cache()


    print(torch.cuda.memory_summary(device = device))

    for i in range(EPOCHS):
        random.shuffle(functions)
        perf = {}
        for func in functions:
            #import pdb; pdb.set_trace()

            print(f"Next dataset: {func['name']}")

            perf[func["name"]] = func["train_func"](i)

            print(torch.cuda.memory_summary(device = device))

            # call the garbage collector and empty the cuda cache mayeb this helps.
            print("GC not called")
            #gc.collect()
            #torch.cuda.empty_cache()


        for m in models_to_train:
            wandb.log({f'{m}_log_sd': wandb.Histogram(models[m].assay_log_sigma.weight.mean(axis = 1).detach().cpu())}, step = i, commit = False)

        perf_ts = [v for p, v in zip(perf.keys(), perf.values()) if p.startswith('tabulaSapiens')]
        perf_rest = [v for p, v in zip(perf.keys(), perf.values()) if not p.startswith('tabulaSapiens')]
        print(f"Perf Tabula Sapiens: {perf_ts}")

        if len(perf_ts) > 0:
            perf_all = perf_rest + [np.mean(perf_ts)]
        else:
            perf_all = perf_rest
        print(perf_all)

        wandb.log({'gene_emb_log_sd': wandb.Histogram(gene_emb.emb_log_sigma.weight.mean(axis = 1).detach().cpu()),
                   "aggregated_performance": np.mean(np.fromiter(perf.values(), dtype = 'float')),
                   "aggregated_performance_mean_ts": np.mean(perf_all),
                   "g_mean_performance_mean_ts": scipy.stats.mstats.gmean(perf_all),
                   "perf_dict": perf,
                  }, 
                  step = i, 
                  commit = True)

        if (i % 50 == 0) & (i > 0):
            gene_emb.get_emb_table().to_csv(
                                            exp_path + f'embedding_{i}_epoch.tsv', 
                                            sep = '\t')
            gene_emb.get_log_sigma_table().to_csv(
                                            exp_path + f'embedding_log_sigma_{i}_epoch.tsv', 
                                            sep = '\t')

            for m_name, m in zip(models, models.values()):
                if m_name != "domainAdapt":
                    m.get_assay_log_sigma_table().to_csv(
                                            exp_path + f'embedding_{m_name}_assay_log_sigma_{i}_epoch.tsv', 
                                            sep = '\t')


        torch.cuda.empty_cache()


    gene_emb.get_emb_table().to_csv(
                                        exp_path + 'embedding.tsv',
                                        sep = '\t')

    gene_emb.get_log_sigma_table().to_csv(
                                            exp_path + f'embedding_log_sigma_{i}_epoch.tsv', 
                                            sep = '\t')

    for m_name, m in zip(models, models.values()):
        if m_name != "domainAdapt":
            m.get_assay_log_sigma_table().to_csv(
                                            exp_path + f'embedding_{m_name}_assay_log_sigma.tsv', 
                                            sep = '\t')


   
    string_average_prec = compute_string_average_precision(gene_emb.get_emb_table()) 
    
    wandb.log({"string_average_prec": string_average_prec})

    print("Model training done!!")
    return string_average_prec



###  Create Optuna study 
study = optuna.create_study(
        study_name=f"{'_'.join(models_to_train_short)}_EPOCHS_{EPOCHS}",
        direction='maximize',
        storage=f"sqlite:///{'_'.join(models_to_train_short)}_EPOCHS_{EPOCHS}",
        sampler = optuna.samplers.TPESampler(multivariate = True),
        load_if_exists=True)

study.enqueue_trial({
        "l_depMap": 0.,
        "l_gtex": 0.,
        "l_perturbSeq": 0.,
        "l_protT5": 0.,
        "l_huMap": 0.,
        "l_tabulaSapiens": 0.,
    }, 
    skip_if_exists = True
)

def get_num_trials_done(study):
    trials_df = study.trials_dataframe()
    if trials_df.empty:
        return 0
    else:
        return len(trials_df.query("state == 'COMPLETE'"))

    
while get_num_trials_done(study) < 50:
    study.optimize(objective)
    




