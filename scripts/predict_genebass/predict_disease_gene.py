

import pandas as pd
import polars as pl

import xgboost as xgb
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.utils.fixes import loguniform

label_table = pl.read_csv(snakemake.input['fold_splits_path'], sep = '\t')
label_table = label_table.filter(pl.col("disease") == snakemake.wildcards['disease'])
label_table = label_table.to_pandas().set_index('gene_id')

emb = pd.read_csv(snakemake.input['embedding_path'], sep = '\t').set_index('gene_id')

dt = label_table.join(emb)

# gene_id not used for prediction as it is in the index.
cols_to_drop = ["target",  "train", "test", "fold", "disease", "min_auprc"]#, "short_name"]


#define the search space for alpha and lambda
param_dist = {
    "reg_lambda": loguniform(1e-2, 1e5), 
    "reg_alpha": loguniform(1e-2, 1e5)
}


results_table = []

for fold in label_table['fold'].unique():
    
    dt_fold = dt.query('fold == @fold')
    X_train = dt_fold.query('train == True').drop(cols_to_drop, axis = 1)
    y_train = dt_fold.query('train == True')['target']
    
    X_test = dt_fold.query('test == True').drop(cols_to_drop, axis = 1)
    y_test = dt_fold.query('test == True')['target']
    
    
    reg = xgb.XGBClassifier(tree_method="hist")
    random_search = RandomizedSearchCV(
        reg, param_distributions=param_dist, n_iter=50, refit = True
    )
    
    mod = random_search.fit(
            X_train, y_train
        )

    pred = mod.predict_proba(
            X_test
        )
    
    fold_results = dt_fold.query('test == True')[cols_to_drop]
    fold_results['pred'] = pred[:,1]
    
    results_table.append(fold_results)

results_table = pd.concat(results_table)
results_table.to_csv(snakemake.output['prediction_results_path'], sep = '\t')
