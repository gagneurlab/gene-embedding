from sklearn.metrics import precision_recall_curve, average_precision_score
import pandas as pd
import polars as pl
import numpy as np

# get paths to all prediction results
prediction_results_path = snakemake.input['prediction_results_path']
# inject dict with study_name to number using input function

# get path for result
eval_table_path = snakemake.output['eval_table_path']

# iterate over paths
prediction_result_table = pl.read_csv(prediction_results_path, sep='\t')


def average_precision_score_polars(x):
    return average_precision_score(x.struct.field("target"), x.struct.field("pred"))


if snakemake.wildcards['stratify_num_pub'] == 'True':
    print("The option to stratify by number of publications was removed.")

else:
    prediction_result_table = prediction_result_table.with_columns([
        pl.struct(["target", "pred"]).apply(average_precision_score_polars).over("fold").alias("average_precision"),
    ])

    eval_table = prediction_result_table.select(["fold", "disease", "min_auprc", "average_precision"]).unique()


eval_table = eval_table.with_column(
    pl.lit(snakemake.wildcards['emb']).alias('emb')
)

eval_table.write_csv(eval_table_path, sep = '\t')







# curve = []
# # iterate over trial per disease
# for name, group in prediction_result_table.groupby(by=['trial', 'disease']):
#     # compute precision-recall
#     prec, recall, _ = precision_recall_curve(group.label, group.prediction)
#     au_prc = average_precision_score(group.label, group.prediction)

#     trial_id = study_root + f"_fold_{select_fold}_{name[0]}"

#     # gather results
#     curve.append(pd.DataFrame({'precision': prec,
#                                'recall': recall,
#                                'fold': select_fold,
#                                'disease': name[1],
#                                'au_prc': au_prc,
#                                'combination': trial_id}))

# # combine curves to one table
# metric_table = pd.concat(curve)
# metric_table['fold'] = metric_table['fold'].astype(str)
# metric_table.to_csv(metric_table_path, index=False, sep='\t')
