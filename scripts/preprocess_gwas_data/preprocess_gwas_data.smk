import pandas as pd
CHROMOSOMES = [str(i) for i in range(1, 23)]

STUDY_TABLE = 'resources/selected_traits.tsv'
#STUDIES = pd.read_csv('scripts/preprocess_gwas_data/Traits_from_Ohler_paper.tsv', sep='\t')['TRAIT'].tolist()

STUDIES = pd.read_csv(STUDY_TABLE, sep='\t')['TRAIT'].tolist()

# STUDIES = ['Triglycerides']

# rule all:
# 	input:
# 		expand("data/magma/{study}/{study}.genes.out", study=STUDIES),
# 		expand("data/magma/{study}/{study}.genes.raw", study=STUDIES),
#         #expand("PoPS/plots/{study}/{study}.{chr}.pops_vs_magma.png", study=STUDIES, chr=CHROMOSOMES)
# 		#expand("PoPS/pops_results/{study}/{study}.{chr}.results", study=STUDIES, chr=CHROMOSOME)
# 		#expand("PoPS/pops_results/{study}.{chr}.coefs", study=STUDIES, chr=CHROMOSOMES)

        

rule get_nealelab_data:
	input:
        	sum_stats_sources = STUDY_TABLE
	params:
        	trait = "{study}"
	output:
        	raw_sum_stats = "data/gwas_data/ukbb_data/{study}/raw_sum_stats.tsv.bgz"
	script:
            "download_sumstats.py"

rule decompress_rawsumstats:
	input:
		raw_sum_stats = "data/gwas_data/ukbb_data/{study}/raw_sum_stats.tsv.bgz"
	output:
		"data/gwas_data/ukbb_data/{study}/raw_sum_stats.tsv"
	shell:
		"""
		bgzip -d {input.raw_sum_stats}
		"""

rule get_snp_annotation:
	params:
		snp_annotation_url = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz"
	output:
		snp_annotation_path_bgzip = "data/gwas_data/ukbb_data/full_variant_qc_metrics.txt.bgz"
	shell:
		"curl {params.snp_annotation_url} -o {output.snp_annotation_path_bgzip}"

rule decompress_snp_annotation:
	input:
		snp_annotation_path_bgzip = "data/gwas_data/ukbb_data/full_variant_qc_metrics.txt.bgz"
	output:
		"data/gwas_data/ukbb_data/full_variant_qc_metrics.txt"
	shell:
		"bgzip -d {input.snp_annotation_path_bgzip}"

rule map_snps:
	input:
		snp_annotation_file = "data/gwas_data/ukbb_data/full_variant_qc_metrics.txt",
		raw_sum_stats = "data/gwas_data/ukbb_data/{study}/raw_sum_stats.tsv",
		#sum_stats_sources = "PoPS/ukbb_data/sum_stats_sources.tsv"
		sum_stats_sources = STUDY_TABLE
	params:
		trait = "{study}"
	output:
		trait_sum_stats = "data/gwas_data/summary_statistics/{study}/sum_stats.tsv"
	threads: 5
	script:
		"map_variants.py"
        
# download reference data as described on the magma homepage: https://ctg.cncr.nl/software/magma        
rule download_reference_data:
    output: 'data/gwas_data/input_data/g1000_eur.zip'
    shell: 'wget -c https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip -O {output}'
     
    
rule extract_reference_data:
    input: 'data/gwas_data/input_data/g1000_eur.zip'
    output: 
        "data/gwas_data/input_data/g1000_eur/g1000_eur.bed",
        "data/gwas_data/input_data/g1000_eur/g1000_eur.bim",
        "data/gwas_data/input_data/g1000_eur/g1000_eur.fam", 
    shell: 'unzip -e {input}'

        
# download gene annotation data as described on the magma homepage: https://ctg.cncr.nl/software/magma
rule download_gene_annot: 
    output: 'data/gwas_data/input_data/NCBI38.zip'
    shell: 'wget -c https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI38.zip -O {output}'
        
rule extract_gene_annotation: 
    input: 'data/gwas_data/input_data/NCBI38.zip'
    output: 'data/gwas_data/input_data/NCBI38/NCBI38.gene.loc'
    shell: "unzip -e {input}"
    
        
			
rule run_magma:
	input:
		"data/gwas_data/input_data/g1000_eur/g1000_eur.bed",
		"data/gwas_data/input_data/g1000_eur/g1000_eur.bim",
		"data/gwas_data/input_data/g1000_eur/g1000_eur.fam",
		gene_annot = "data/gwas_data/input_data/magma_0kb.genes.annot",
#		gene_annot = 'data/gwas_data/input_data/NCBI38/NCBI38.gene.loc',
		pval = "data/gwas_data/summary_statistics/{study}/sum_stats.tsv",
#		meta_file = "PoPS/summary_statistics/{study}/meta_information.txt"
	output:
		"data/magma/{study}/{study}.genes.out",
		"data/magma/{study}/{study}.genes.raw"
	threads: 10
	shell:
		"magma_v1.10_static/magma "
		"--bfile data/gwas_data/input_data/g1000_eur/g1000_eur "
		"--gene-annot {input.gene_annot} "
		"--pval {input.pval} ncol=N "
		"--gene-model snp-wise=mean "
		"--out data/magma/{wildcards.study}/{wildcards.study}"
		

## before running POPS one needs to download the data from here: 
# https://www.finucanelab.org/data
# https://www.dropbox.com/sh/o6t5jprvxb8b500/AADZ8qD6Rpz4uvCk0b5nUnPaa/data?dl=0&subfolder_nav_tracking=1

rule run_pops_feature_selection:
	input:
		"data/magma/{study}/{study}.genes.out",
		"data/magma/{study}/{study}.genes.raw",
		feature = "data/input_data/PoPS/PoPS.features.txt.gz"
	output:
		"data/PoPS/pops_results/{study}/{study}.features"
	threads: 5
	shell:
		"python scripts/preprocess_gwas_data/pops.feature_selection.py "
		"--features {input.feature} "
		"--gene_results data/magma/{wildcards.study}/{wildcards.study} "
		"--out data/PoPS/pops_results/{wildcards.study}/{wildcards.study}"

rule run_pops_prediction:
	input:
		"data/magma/{study}/{study}.genes.out",
		"data/magma/{study}/{study}.genes.raw",
		selected_features = "data/PoPS/pops_results/{study}/{study}.features",
		features = "data/input_data/PoPS/PoPS.features.txt.gz",
		gene_loc = "data/input_data/PoPS/gene_loc.txt",
		control_features = "data/input_data/PoPS/control.features",
	output:
		"data/PoPS/pops_results/{study}/{study}.{chr}.results",
		"data/PoPS/pops_results/{study}/{study}.{chr}.coefs"
	threads: 5
	shell:
		"python scripts/preprocess_gwas_data/pops.predict_scores.py "
		"--gene_loc {input.gene_loc} "
		"--gene_results data/magma/{wildcards.study}/{wildcards.study} "
		"--features {input.features} "
		"--selected_features {input.selected_features} "
		"--control_features {input.control_features} "
		"--chromosome {wildcards.chr} "
		"--out data/PoPS/pops_results/{wildcards.study}/{wildcards.study}"
