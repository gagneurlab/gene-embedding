

import itertools

list_individual_embeddings.sort()

all_combination = list(itertools.chain.from_iterable(
    [itertools.combinations(list_individual_embeddings, r = i) for i in range(1, len(list_individual_embeddings) + 1)]
))
all_combination = {'.'.join(i) : i for i in all_combination}

all_combination_list = [f"combined_{e}" for e in all_combination.keys()]



rule create_combinations:
    params:
        norm = "T"
    input: embedding_paths = lambda wildcards: [config['datadir'] + '/embedding/individual/' + e + '.tsv' for e in all_combination[wildcards.emb]]
    output: out_path = config['datadir'] + '/embedding/combination/combined_{emb}.tsv'
    script: "create_combinations.py"
        
    

