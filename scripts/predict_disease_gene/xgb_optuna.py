#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import precision_recall_curve, average_precision_score
import sklearn.model_selection
import optuna
from scipy import stats
import os

np.random.seed(1234)


# ## Explore Hyper Params With Optuna
# change this parameter for new study
trial_name = snakemake.wildcards.study_root


# ### Define Model Loop
def compute(embedding_df, coding_gene_table, fold_splits, max_iter, max_depth, learning_rate, verbose=False):
    curve = []
    prediction_table = []
    # things that do not need to be done for each disease
    
#     # make sure we are only using protein coding genes
#     embedding_df = embedding_df.merge(coding_gene_table[['Gene stable ID']], left_on='gene_id', right_on='Gene stable ID').drop('Gene stable ID', axis=1).set_index('gene_id')
    
#     # get difference to universe
#     gene_count_difference = len(coding_gene_table)-len(embedding_df)

    # make sure we are only using protein coding genes
    embedding_df = embedding_df[embedding_df.index.isin(coding_gene_table['Gene stable ID'])]
    
    # get difference to universe
    gene_count_difference = len(coding_gene_table['Gene stable ID'].unique())-len(embedding_df)
    
    vali_curve = []
    test_predictions = []
    
    # group by all disease and iterate over them
    # there will be only one fold here at a time buut code is okay as is
    for group_tuple, disease_group in fold_splits.groupby(by=['disease', 'fold']):
        disease_name = group_tuple[0]
        fold = group_tuple[1]

        # define model
        boost = HistGradientBoostingClassifier(max_iter=max_iter, max_depth=max_depth, learning_rate=learning_rate)

        # merge for embedding for selected trait
        combined_table = embedding_df.join(disease_group, how='inner')
        
        # query for true test and true target label and sample 50%
        test_true_combined = combined_table.query('test == True and target == True').sample(frac=0.5, random_state=1234)
        # query for true test and false target label
        test_false_combined = combined_table.query('test == True and target == False').sample(frac=0.5, random_state=1234)
        # combine new test split
        test_combined_table = pd.concat([test_true_combined, test_false_combined])
        
        # create new test split by droping val split genes
        val_combined_table = combined_table.drop(test_combined_table.index)
        
        if verbose:
            print(disease_name + f',fold {fold}')
            print(len(val_combined_table))
            print(len(test_combined_table))

        # subset the embedding df
        train_emb = val_combined_table.query('train')[embedding_df.columns]

        # train random forest                           
        boost.fit(train_emb.values, val_combined_table.query('train').target)

        val_emb = val_combined_table.query('test')[embedding_df.columns]
        # compute prediction and get target class
        val_pred = boost.predict_proba(val_emb)[:,1]

        # fill up all arrays to be on the same universe
        label_array = np.append(val_combined_table.query('test').target, np.repeat(False, gene_count_difference))
        prediction_array = np.append(val_pred, np.repeat(0, gene_count_difference))
        gene_array = np.append(val_combined_table.query('test').index, np.repeat('dummy', gene_count_difference))

        # compute precision-recall
        #prec, recall, _ = precision_recall_curve(val_combined_table.query('test').target, val_pred)
        au_prc = average_precision_score(val_combined_table.query('test').target, val_pred)

        # gather precision-recall results
        curve.append(pd.DataFrame({'fold': fold, 
                                   'disease': disease_name,
                                   'au_prc': au_prc}, index=[0]))

        # run predictions on new test set
        test_emb = test_combined_table.query('test')[embedding_df.columns]
        # compute prediction and get target class
        test_pred = boost.predict_proba(test_emb)[:,1]

        # fill up all arrays to be on the same universe
        label_array = np.append(test_combined_table.query('test').target, np.repeat(False, gene_count_difference))
        prediction_array = np.append(test_pred, np.repeat(0, gene_count_difference))
        gene_array = np.append(test_combined_table.query('test').index, np.repeat('dummy', gene_count_difference))

        # gather prediction results
        prediction_table.append(pd.DataFrame({
            'gene_id': gene_array,
            'fold': fold,
            'prediction': prediction_array,
            'label': label_array,
            'disease': disease_name,
            'embedding': trial_name
        }))


    # combine results
    # precision-recall
    prcurve = pd.concat(curve)
    prcurve['fold'] = prcurve['fold'].astype(str)
    
    # combine predictions from all disaeses
    prediction_table = pd.concat(prediction_table)

    # make fold string to be categorical
    prediction_table['fold'] = prediction_table['fold'].astype(str)

    return prcurve, prediction_table

# (https://optuna.readthedocs.io/en/stable/faq.html#objective-func-additional-args).
def objective(trial):
    print(f'Trial Number {trial.number}')
    # recomended to tune for xdg:
    #Number of trees.
    #Tree depth.
    #Step Size (learning rate).
    learning_rate = trial.suggest_float('learnign_rate', 0.001, 1, log=True)
    max_iter = trial.suggest_int('max_iter', 10, 10 * len(embedding_df.columns), log=True)
    max_depth = trial.suggest_int('max_depth', 1, 10, log=True)
    
    val_prcurve, test_prediction_table = compute(embedding_df, coding_gene_table, fold,
                                            max_iter=max_iter,
                                            max_depth=max_depth,
                                            learning_rate=learning_rate
                                           )
    
    test_prediction_table['trial'] = trial.number
    test_prediction_table['study'] = per_fold_study
    
    # save prediction results to file
    output_path = test_predictions_folder + trial_name + '.tsv'
    test_prediction_table.to_csv(output_path, sep='\t', index=False, header=not os.path.exists(output_path) ,mode='a')
    
    return val_prcurve.groupby('disease').au_prc.median().median()

def get_num_trials_done(study):
    trials_df = study.trials_dataframe()
    if trials_df.empty:
        return 0
    else:
        return len(trials_df.query("state == 'COMPLETE'"))

# ### Load Data
# keep this in case this becomes part of snakemake pipeline
#embedding_path = '../data/AE_Embeddings/combinations/pca-gtex_pca-crispr_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall.tsv'
#embedding_path = '../data/AE_Embeddings/combinations/gtex_crisprembnew_stringnolitverse_proteinsmall.tsv'
disease_gene_path = 'data/input_data/disease_gene_eval/disease_gene_table.tsv'
coding_gene_path = 'data/input_data/disease_gene_eval/ensembl_id_map.tsv'
fold_splits_path = "data/modeling_splits/disease_gene_splits.tsv"
test_predictions_folder = 'data/Evaluation/hyperopt/predictions/'

# load main data 
embedding_df = pd.read_csv(snakemake.input['embedding_path'], sep='\t').set_index('gene_id')
disease_gene_table = pd.read_csv(disease_gene_path, sep='\t')
coding_gene_table = pd.read_csv(coding_gene_path, sep='\t').dropna()
fold_splits = pd.read_csv(fold_splits_path, sep='\t', index_col='gene_id')


# select traits with sufficient genes
df = disease_gene_table.disease.value_counts().rename_axis('disease').reset_index(name='counts')
disease_gene_table = disease_gene_table.merge(df, on='disease')
larger_traits = disease_gene_table.query('counts >= 50').disease.unique()
print(larger_traits)
fold_splits = fold_splits.query('disease in @larger_traits')


#get name of embedding from path

for f in fold_splits.fold.unique():

    fold = fold_splits.query('fold == @f')
    
    per_fold_study = trial_name + f'_fold_{f}'
    
    study = optuna.create_study(
            study_name= per_fold_study,
            direction = 'maximize',
            storage=f"sqlite:///xgb_opt_{snakemake.wildcards['study_root'].split('.')[0]}.db",
            sampler = optuna.samplers.TPESampler(multivariate = True),
            load_if_exists=True)

    while (get_num_trials_done(study) < 100):
        study.optimize(objective, n_trials=1)
