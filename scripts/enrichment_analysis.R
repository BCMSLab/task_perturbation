# Perform enrichment analysis
#
# Inputs: hallmarks gene sets, archetypes of k=6, and drug response
# Outputs: archetype enrichment and response enrichment

# load libraries
library(tidyverse)
library(reshape2)
library(archetypes)
library(SummarizedExperiment)
library(clusterProfiler)

# enrichment, by archetype
hallmarks <- read_csv('output/hallmarks_gene_sets.csv')
names(hallmarks) <- c('gene', 'term')
hallmarks <- dplyr::select(hallmarks, term, gene)

a6 <- read_rds('output/archetypes_6k.rds')
arcs <- t(parameters(a6))
colnames(arcs) <- 1:ncol(arcs)

arc_list <- lapply(seq_len(ncol(arcs)), function(x) arcs[, x])

stats <- map(arc_list, function(x) {
  vec <- scale(x)
  names(vec) <- rownames(arcs)
  sort(vec, decreasing = TRUE)
})

set.seed(123)
enrich_archetype <- map_df(stats, function(x) {
  enrich <- GSEA(x,
                 pAdjustMethod = 'fdr',
                 TERM2GENE = as.data.frame(hallmarks),
                 pvalueCutoff = 1,
                 maxGSSize = 1500,
                 nPermSimple = 10000)
  as_tibble(enrich)
}, .id = 'archetype')

write_csv(enrich_archetype, 'output/archetype_enrichment.csv')

# # enrichment, by arch in task
# relevant_terms <- read_csv('data_clean/relevant_terms.csv') %>% na.omit()
# task_gene <- left_join(term_gene, relevant_terms, by = c('term'='ID')) %>%
#   na.omit() %>%
#   dplyr::select(term = task, gene) %>%
#   unique()
# 
# set.seed(123)
# enrich_archetype_task <- map_df(stats, function(x) {
#   enrich <- GSEA(x,
#                  pAdjustMethod = 'fdr',
#                  TERM2GENE = task_gene,
#                  pvalueCutoff = 1)
#   as_tibble(enrich)
# }, .id = 'archetype')
# 
# write_rds(enrich_archetype_task, 'data_clean/enrich_archetype_task.rds')

# enrichment, by treatment
drug_response <- read_rds('output/drug_response.rds')

response_stats <- split(assay(drug_response),
                        rep(1:ncol(drug_response),
                            each = nrow(drug_response)))
response_stats <- map(response_stats, function(x) {
  names(x) <- rownames(drug_response)
  sort(x, decreasing = TRUE)
})

names(response_stats) <- colnames(drug_response)

stats <- map(stats, function(x) x[!duplicated(names(x))])

set.seed(123)

enrich_treatment <- map_df(response_stats, function(x) {
  enrich <- GSEA(x,
                 pAdjustMethod = 'fdr',
                 TERM2GENE = hallmarks,
                 pvalueCutoff = 1,
                 maxGSSize = Inf)
  as_tibble(enrich)
}, .id = 'group') %>%
  separate(group, c('cell_id', 'drug'))

write_csv(enrich_treatment, 'output/response_enrichment.csv')

# # enrichment, by treatment in task
# set.seed(123)
# enrich_treatment_task <- map_df(stats, function(x) {
#   enrich <- GSEA(x,
#                  pAdjustMethod = 'fdr',
#                  TERM2GENE = task_gene,
#                  pvalueCutoff = 1)
#   as_tibble(enrich)
# }, .id = 'group') %>%
#   separate(group, c('cell_id', 'drug'))
# 
# write_rds(enrich_treatment_task, 'data_clean/enrich_treatment_task.rds')
