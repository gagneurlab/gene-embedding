import pandas as pd

per_fold_results_path = snakemake.input['per_fold_results_path']
best_prerec_path = snakemake.output['best_table_path']
study_root = snakemake.wildcards['study_root']

study_result_table = pd.concat([pd.read_csv(p, sep='\t') for p in per_fold_results_path])

study_result_table.groupby(by='combination').au_prc.median().max()

best_folds_table = pd.DataFrame(study_result_table.groupby('combination').au_prc.median()).reset_index().sort_values(by='au_prc', ascending=False)

best_folds_table['study'] = best_folds_table.combination.str.rsplit(pat="_" ,n=1)
best_folds_table['study'] = best_folds_table.study.apply(lambda x: x[0])

best_folds_table = best_folds_table.drop_duplicates(subset=['study'])

tmp_table = []
for c in best_folds_table.combination:
    print(c)
    tmp_table.append(study_result_table.query('combination == @c'))

best_prerec_table = pd.concat(tmp_table)
best_prerec_table['trial_number'] = best_prerec_table.combination.str.rsplit(pat="_" ,n=1)
best_prerec_table['trial_number'] = best_prerec_table.trial_number.apply(lambda x: x[1])
best_prerec_table['combination'] = study_root

best_prerec_table.to_csv(best_prerec_path, index=False, sep='\t')
