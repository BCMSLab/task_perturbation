#!/bin/bash

# This file runs inside the bcmslab/task_perturbation2

all: output/hallmarks_gene_sets.csv \
		output/growth_rates.csv

# cancer hallmarks
output/hallmarks_gene_sets.csv: scripts/hallmarks_gene_sets.R
		R CMD BATCH --vanilla scripts/hallmarks_gene_sets.R

# growth rates
output/growth_rates.csv: scripts/growth_rates.R 
		R CMD BATCH --vanilla scripts/growth_rates.R

# gene expression
output/gene_expression.rds: scripts/gene_expression.R \
		output/growth_rates.csv
		R CMD BATCH --vanilla scripts/gene_expression.R

# drug response (fold-changes)
output/drug_response.rds: scripts/drug_response.R \
		output/gene_expression.rds
		R CMD BATCH --vanilla scripts/drug_response.R
		
# archetype analysis
output/archtypes.rds: scripts/archetype_analysis.R \
		output/gene_expression.rds
		R CMD BATCH --vanilla scripts/archetype_analysis.R

# variance partioning
output/variance_partition.csv: scripts/variance_partition.R \
		output/gene_expression.rds \
		output/gene_expression_averages.csv \
		output/distances.csv
		R CMD BATCH --vanilla scripts/variance_partition.R

# multi dimension scaling
output/dimension_scaling.csv: scripts/dimension_scaling.R \
		output/gene_expression.rds \
		output/gene_expression_averages.csv \
		output/distances.csv
		R CMD BATCH --vanilla scripts/dimension_scaling.R
		
# enrichment analysis
output/archetype_enrichment.csv: scripts/enrichment_analysis.R \
		output/hallmarks_gene_sets.csv \
		output/archetypes_6k.rds \
		output/drug_response.rds
		R CMD BATCH --vanilla scripts/enrichment_analysis.R
		
# contextual analysis
output/mod_4.csv: scripts/contextual_analysis.R \
		output/hallmarks_gene_sets.csv \
		output/drug_response.rds \
		output/gene_expression_averages.csv
		R CMD BATCH --vanilla scripts/contextual_analysis.R
