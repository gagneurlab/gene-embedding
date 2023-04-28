import polars as pl
import pandas as pd
import numpy as np
# import xgboost as xgb
# import plotnine as p9
import scipy.stats
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, GroupKFold

from sklearn.linear_model import ElasticNetCV
from sklearn.datasets import make_regression

from sklearn.utils.fixes import loguniform

# adapded from pops 
def build_control_covariates(metadata):
	genesize = metadata.NPARAM.values.astype(float)
	genedensity = metadata.NPARAM.values/metadata.NSNPS.values
	inverse_mac = 1.0/metadata.MAC.values
	cov = np.stack((genesize, np.log(genesize), genedensity, np.log(genedensity), inverse_mac, np.log(inverse_mac)), axis=1)
	return cov

def munge_sigma(magma_gene_raw):
	f = open(magma_gene_raw)
	lines = list(f)[2:]
	lines = [np.asarray(line.strip('\n').split(' ')) for line in lines]
	sigmas = []
	gene_metadata = []
	gene_lists = []
	for chrom in range(1,23):
		chr_start = min(np.where([int(line[1])==chrom for line in lines])[0])
		chr_end = max(np.where([int(line[1])==chrom for line in lines])[0])
		lines_chr = lines[chr_start:chr_end+1]
		n_genes = len(lines_chr)
		sigma_chr = np.zeros([n_genes, n_genes])
		gene_NSNPs = np.zeros(n_genes)
		gene_NPARAM = np.zeros(n_genes)
		gene_MAC = np.zeros(n_genes)
		for i in range(n_genes):
			line = lines_chr[i]
			gene_NSNPs[i] = line[4]
			gene_NPARAM[i] = line[5]
			gene_MAC[i] = line[7]
			if line.shape[0] > 9:
				gene_corrs = np.asarray([float(c) for c in line[9:]])
				sigma_chr[i, i-gene_corrs.shape[0]:i] = gene_corrs
		sigma_chr = sigma_chr+sigma_chr.T+np.identity(n_genes)
		sigmas.append(sigma_chr)
		gene_metadata_chr = pd.DataFrame(data={'NSNPS': gene_NSNPs, 'NPARAM': gene_NPARAM, 'MAC': gene_MAC})
		gene_metadata.append(gene_metadata_chr)
		gene_list_chr = [line[0] for line in lines_chr]
		gene_lists.append(gene_list_chr)
	return sigmas, gene_metadata, gene_lists


def compute_Ls(sigmas, Y):
	Ls = []
	min_lambda=0
	for sigma in sigmas:
		W = np.linalg.eigvalsh(sigma)
		min_lambda = min(min_lambda, min(W))
	#Y = pd.read_table(args.gene_results+'.genes.out', delim_whitespace=True).ZSTAT.values
	ridge = abs(min_lambda)+.05+.9*max(0, np.var(Y)-1)
	for sigma in sigmas:
		sigma = sigma+ridge*np.identity(sigma.shape[0])
		L = np.linalg.cholesky(np.linalg.inv(sigma))
		Ls.append(L)
	return Ls


# sigmas, metadata, gene_lists = munge_sigma('../data/PoPS/magma_results/RBC/RBC.genes.raw')
sigmas, metadata, gene_lists = munge_sigma(snakemake.input.magma_gene_raw)

chrom = int(snakemake.wildcards['chrom'])


# create covariates from pops
covariates = []
for i in range(0, 22):
    print(i)
    covariates.append(pd.DataFrame(build_control_covariates(metadata[i]),
                                   index = gene_lists[i], 
                                   columns = ['genesize',
                                              'log_genesize',
                                              'genedensity',
                                              'log_genedensity',
                                              'inverse_mac',
                                              'log_inverse_mac'])
                      )


covariates = pd.concat(covariates)


# load embedding
emb = pd.read_csv(snakemake.input['embedding_path'], sep = '\t').set_index('gene_id')
# emb = pd.read_csv(snakemake.input.emb, sep = "\t").set_index("ENSGID")
#emb = pd.read_csv('../data/AE_Embeddings/combinations/norm-gtex_norm-crisprembnew_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall.tsv', 
#                  sep = "\t").set_index("gene_id")
emb

magma = pd.read_csv(snakemake.input.magma_output, delim_whitespace=True)
# magma = pd.read_csv('../data/PoPS/magma_results/RBC/RBC.genes.out', delim_whitespace=True)

magma = magma.merge(covariates, left_on = "GENE", right_index = True)



# project Y to LY
## compute Ls
Ls = compute_Ls(sigmas, magma.ZSTAT)

def project_Y(Ls, magma_Z):
    LYs = []
    for i in range(22):
        L = Ls[i]
        magma_temp = magma.set_index("GENE").reindex(gene_lists[i]).reset_index()
        
        LYs.append(pd.DataFrame({"GENE": magma_temp.GENE, "LY": np.matmul(L, magma_temp.ZSTAT)}))
    return pd.concat(LYs)

def project_Y_back(Ls, res):
    LYs = []
    for i in range(22):
        L = np.linalg.inv(Ls[i])
        temp = res.set_index("GENE").reindex(gene_lists[i]).reset_index()
        
        LYs.append(pd.DataFrame({"GENE": temp.dropna().GENE,
                                 "pred": np.matmul(L[~temp.pred_LY.isna(), :][:, ~temp.pred_LY.isna()], 
                                                   temp.dropna().pred_LY), 
                                }))
    return pd.concat(LYs)


def project_Y_back_chrom(Ls, res, chrom):
        L = np.linalg.inv(Ls[int(chrom) - 1])
        temp = res.set_index("GENE").reindex(gene_lists[int(chrom) -1]).reset_index()
        
        LYs = pd.DataFrame({"GENE": temp.dropna().GENE,
                                 "pred": np.matmul(L[~temp.pred_LY.isna(), :][:, ~temp.pred_LY.isna()], 
                                                   temp.dropna().pred_LY), 
                                })
        return LYs



magma = magma.merge(project_Y(Ls, magma))

# merge with embedding
dt = magma.merge(emb, left_on = "GENE", right_on = "gene_id")
# dt = magma.merge(emb, left_on = "GENE", right_on = "ENSGID")

# reg_lambda = 1000, reg_alpha = 100
param_dist = {
    "reg_lambda": loguniform(1e-2, 1e5),
    "reg_alpha": loguniform(1e-2, 1e5)
}

n_iter_search = 20

# Utility function to report best scores
# def report(results, n_top=3):
#     for i in range(1, n_top + 1):
#         candidates = np.flatnonzero(results["rank_test_score"] == i)
#         for candidate in candidates:
#             print("Model with rank: {0}".format(i))
#             print(
#                 "Mean validation score: {0:.3f} (std: {1:.3f})".format(
#                     results["mean_test_score"][candidate],
#                     results["std_test_score"][candidate],
#                 )
#             )
#             print("Parameters: {0}".format(results["params"][candidate]))
#             print("")


# ----------------- XGBoost -----------------
# reg = xgb.XGBRegressor(tree_method="hist")
# gkf = GroupKFold(n_splits=5)
# random_search = RandomizedSearchCV(
#     reg, param_distributions=param_dist, n_iter=n_iter_search, refit = True, 
#     cv = gkf.split(dt.query("CHR != @chrom")['LY'], 
#                    dt.query("CHR != @chrom")['LY'], 
#                    groups = dt.query("CHR != @chrom")['CHR']
#                    )
# )
# mod = random_search.fit(
#     dt.query("CHR != @chrom").drop(["GENE", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N", "ZSTAT", "P", "LY"], axis=1),
#     dt.query("CHR != @chrom")['LY']
# )

# ----------------- Elastic Net -----------------
gkf = GroupKFold(n_splits=5)
reg = ElasticNetCV(random_state=0, 
                   cv=gkf.split(dt.query("CHR != @chrom")['LY'], 
                                dt.query("CHR != @chrom")['LY'], 
                                groups = dt.query("CHR != @chrom")['CHR']
                               )
                  )

mod = reg.fit(dt.query("CHR != @chrom").drop(["GENE", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N", "ZSTAT", "P", "LY"], axis=1),
              dt.query("CHR != @chrom")['LY']
             )

pred = mod.predict(
    dt.query("CHR == @chrom").drop(["GENE", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N", "ZSTAT", "P", "LY"], axis=1),
)

df_chrom = dt.query("CHR == @chrom")[["GENE", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N", "ZSTAT", "P", "LY"]]
df_chrom['pred_LY'] = pred


# print(f"Chrom: {chrom}: R2: {mod.score(df_chrom.LY, df_chrom.pred_LY)[0]**2}")
print(f"Chrom: {chrom}: R2: {scipy.stats.pearsonr(df_chrom.LY, df_chrom.pred_LY)[0]**2}")
print("------------- CV results ------------------")
# print(report(random_search.cv_results_))
print("------------- CV results END --------------")

print(df_chrom)


print(project_Y_back_chrom(Ls, df_chrom, chrom))


df_chrom = df_chrom.merge(project_Y_back_chrom(Ls, df_chrom, chrom))


# print(f"Overall R2: {scipy.stats.pearsonr(df.ZSTAT, df.pred)[0]**2}")

# print("Per chrom R2: \n")

# for i in range(1, 23):
#     df_tmp = df.query("CHR == @i")
#     print(scipy.stats.pearsonr(df_tmp.ZSTAT, df_tmp.pred)[0]**2)

df_chrom.to_csv(snakemake.output.pred, sep = '\t')
