# Perform archetype analysis
#
# Inputs: gene expression
# Outputs: archetypes, archetypes of k=6, and distances 

# load libraries
library(tidyverse)
library(reshape2)
library(archetypes)
library(SummarizedExperiment)
library(clusterProfiler)

# load cancer cell genes
gene_expression <- read_rds('output/gene_expression.rds')
tissues <- c('breast', 'large_intestine', 'ovary', 'lung',
             'haematopoietic_and_lymphoid_tissue')
se <- gene_expression[, gene_expression$tissue %in% tissues]

dat <- aggregate(t(assay(se)),
                 by = list(cell_id = se$cell_id),
                 FUN = mean)

dat <- column_to_rownames(dat, 'cell_id')
dat <- dat[, !is.na(names(dat))]

write.csv(as.data.frame(dat), 'output/gene_expression_averages.csv')

## perform archetype analysis
set.seed(123)
as <- stepArchetypes(dat,
                     k = 1:10,
                     verbose = TRUE)

write_rds(as, 'output/archetypes.rds')

a6 <- bestModel(as[[6]])

write_rds(a6, 'output/archetypes_6k.rds')

# distances
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

distances <- melt(distances)
names(distances) <- c('cell_id', 'archetype', 'dss')

write_csv(distances, 'output/cell_distances.csv')
