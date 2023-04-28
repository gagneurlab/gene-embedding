# A genome-wide experiment-based functional gene embedding

The notion of gene function is central to genetics. However, functional annotations that are free text or based on controlled vocabularies are complex and often partial or missing. Therefore, their integration into machine learning models is cumbersome and biassed for well-studied genes.

Functional gene embeddings, numerical vectors capturing gene function, offer an exciting avenue to address this issue. Such embeddings are obtained by applying self-supervised approaches on various data types including quantitative measurements, protein interaction networks, and literature. However, their utility and relative performance for diverse prediction tasks has not been assessed so far. 


This repo contains pipelines to create all the embeddings and figures described in

>>>> add link here <<<<

# Cloning the pipeline 

To run the pipeline on your sit first clone the repository to your system.
The repository further depends on node2vec and VERSE. Checkout the two submodules using:
```
git submodule update --init
```
After cloning the submodules the need to get compiled. Depeding on where they are compiled the paths in the sankemake file have to get adjusted.


# The analysis scripts are structured as follows

The scripts folders contain folders for each step of the pipeline. 
The full pipeline can be triggered by running the global Snakefile using 

```
snakemake -j 5 ...
```

The pipeline is structurred into preprocessing scripts, embedding scripts and prediction scripts.

The prediction scripts are split into folders depending on the tasks they are trained on.

## Disease gene prediction

For six curated disease-gene lists train an XGBoost model in five-fold cross-validation, predicting for each gene whether it is part of each list. For each disease, all non-annotated genes are treated as negative instances.


## HPO prediction

A simple two-layer multiclass neural network to predict for each gene it's association with any of the selected [HPO](https://hpo.jax.org/app/) terms. HPO terms associated with at least genes are selected, resulting a set of 1994 terms. Thirty percent of the genes are used to evaluate the performance of the model.

## Cancer prediction
Association of genes with cancer are obtained from [TCGA](https://www.cancer.gov/ccg/research/genome-sequencing/tcga). An XGBoost model learns to predict if a gene is associated with cancer and the model is evaluated using five-fold cross cross-validation.


## RVAS prediction

Predict trait-gene associations reported by [genebass](https://app.genebass.org/), a rare variant association study performed on  the UK Biobank. Traits associated with more than 30 genes are selected. Further filtering is performed to remove traits which have disproportionate number of overlapping genes, resulting in a set of 34 genes. An XGBoost model learns the associations and five-fold cross-validation is used for evaluating the performance.


## GWAS prediction

Predict gene-level association scores for 22 uncorrelated blood biomarker GWAS studies from the [Pan-UK Biobank](https://pan.ukbb.broadinstitute.org/) study. For a given trait, [MAGMA](https://ctg.cncr.nl/software/magma) summarises gene-trait associations in one z-score per gene. An Elastic Net model learns to predict the projected z-score from the functional gene embeddings and six additional covariates. Leave-One-Chromosome-Out strategy used for evaluating the performance.




