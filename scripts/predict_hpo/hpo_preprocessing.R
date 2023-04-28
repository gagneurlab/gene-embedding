

library(data.table)
library(tidyr)


# code used to format the header fread friendly.

# gsub('<tab>', ' ', 
#      gsub(' ', '_',
#           gsub('-', '_', "entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-ID<tab>HPO-Term-Name<tab>Frequency-Raw<tab>Frequency-HPO<tab>Additional Info from G-D source<tab>G-D source<tab>disease-ID for link")
#      )
# )

print(snakemake@input['gene2pheno'])
hpo_gene <- fread(snakemake@input[['gene2pheno']])

# Number of unique HPO terms
hpo_gene[, uniqueN(HPO_Term_ID)]

# Number of annotated genes
hpo_gene[, uniqueN(entrez_gene_id)]

# Numer of genes per HPO term
hpo_gene[, uniqueN(entrez_gene_id), by = HPO_Term_ID][, V1] %>% ecdf() %>% plot()


hpo_gene[, uniqueN(entrez_gene_id), by = HPO_Term_ID][order(V1)]

# There are 1995 HPO terms with at least 20 genes annotated.
hpo_gene[, uniqueN(entrez_gene_id), by = HPO_Term_ID][order(V1)][V1 >= 20]


gene_mapping <- fread(snakemake@input[['gene_mapping']])
gene_mapping
colnames(gene_mapping) <- gsub(' ', '_', colnames(gene_mapping))
gene_mapping

#gene_mapping <- separate(gene_mapping, HGNC_ID, into = c('hgnc', 'hgnc_id'))

gene_mapping <- unique(gene_mapping[, .(Approved_symbol,
                                        Ensembl_ID = `Ensembl_ID(supplied_by_Ensembl)`,
                                        NCBI_Gene_ID = `NCBI_Gene_ID(supplied_by_NCBI)`
                                        )])
#gene_mapping[, hgnc_id := as.integer(hgnc_id)]
gene_mapping

hpo_gene

hpo_gene_ensmbl <- merge(hpo_gene, gene_mapping, by.x = 'entrez_gene_id', by.y = 'NCBI_Gene_ID')

hpo_gene_ensmbl[, uniqueN(Ensembl_ID), by = HPO_Term_ID][, V1 ] %>% ecdf() %>% plot()
hpo_gene_ensmbl[, uniqueN(Ensembl_ID), by = HPO_Term_ID][  V1 >= 20 ] 

hpo_gene_ensmbl[, n_genes_annotated := uniqueN(Ensembl_ID), by = HPO_Term_ID]
hpo_gene_ensmbl[, .(N = uniqueN(Ensembl_ID)), by = Ensembl_ID][, N]  %>% ecdf() %>% plot()

hpo_gene_ensmbl[, n_hpo_terms_per_gene := uniqueN(HPO_Term_ID), by = Ensembl_ID]

fwrite(hpo_gene_ensmbl, snakemake@output[['gene2pheno']], sep = '\t')

