# load libraries
library(tidyverse)
library(reshape2)
library(archetypes)
library(SummarizedExperiment)
library(clusterProfiler)

# load cancer cell genes
se <- read_rds('data_clean/cancer_cells_genes.rds')
tissues <- c('breast', 'large_intestine', 'ovary', 'lung',
             'haematopoietic_and_lymphoid_tissue')
se <- se[, se$tissue %in% tissues]

dat <- aggregate(t(assay(se)),
                 by = list(cell_id = se$cell_id),
                 FUN = mean)

dat <- column_to_rownames(dat, 'cell_id')
dat <- dat[, !is.na(names(dat))]

write_rds(dat, 'data_clean/cancer_cells_genes_average.rds')

## perform archetype analysis
set.seed(123)
as <- stepArchetypes(dat,
                     k = 1:10,
                     verbose = TRUE)
screeplot(as)
write_rds(as, 'data_clean/cell_archetypes.rds')

a6 <- bestModel(as[[6]])
simplexplot(a6)

write_rds(a6, 'data_clean/archetypes_6.rds')

# distances
a6 <- read_rds('data_clean/archetypes_6.rds')
arcs <- t(parameters(a6))
colnames(arcs) <- 1:ncol(arcs)

arc_list <- lapply(seq_len(ncol(arcs)), function(x) arcs[, x])

dat_mat <- as.matrix(dat)

distances <- lapply(arc_list,
                    function(a) {
                      diff <- a - t(dat_mat)
                      sqrt(colSums(diff^2))
                    })

distances <- as.matrix(bind_cols(distances))
colnames(distances) <- 1:ncol(distances)
rownames(distances) <- rownames(dat)
write_rds(distances, 'data_clean/arc_cell_distances.rds')

# enrichment, by archetype
term_gene <- read_rds('data_clean/hallmarks_sets.rds') %>%
  melt() %>%
  setNames(c('gene', 'term')) %>%
  dplyr::select(term, gene)

stats <- map(arc_list, function(x) {
  vec <- scale(x)
  names(vec) <- rownames(arcs)
  sort(vec, decreasing = TRUE)
})

set.seed(123)
enrich_archetype <- map_df(stats, function(x) {
  enrich <- GSEA(x,
                 pAdjustMethod = 'fdr',
                 TERM2GENE = term_gene,
                 pvalueCutoff = 1,
                 maxGSSize = 1500,
                 nPermSimple = 10000)
  as_tibble(enrich)
}, .id = 'archetype')

write_rds(enrich_archetype, 'data_clean/enrich_archetype.rds')

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
cancer_cells_scores <- read_rds('data_clean/cancer_cells_scores.rds')
stats <- cancer_cells_scores %>%
  imap(function(x, .y) {
    drug <- str_split(.y, '_', simplify = TRUE)[, 2]
    vec <- unlist(x[, drug])
    names(vec) <- x$symbol
    sort(vec, decreasing = TRUE)
  })
stats <- map(stats, function(x) x[!duplicated(names(x))])

set.seed(123)
enrich_treatment <- map_df(stats, function(x) {
  enrich <- GSEA(x,
                 pAdjustMethod = 'fdr',
                 TERM2GENE = term_gene,
                 pvalueCutoff = 1,
                 maxGSSize = Inf)
  as_tibble(enrich)
}, .id = 'group') %>%
  separate(group, c('cell_id', 'drug'))

write_rds(enrich_treatment, 'data_clean/enrich_treatment.rds')

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
