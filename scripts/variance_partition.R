# Partition the variance in the typical gene expression
#
# Inputs: gene expression averages, archetype distances
# Output: var part object

# load libraries
library(tidyverse)
library(variancePartition)

# Intialize parallel backend with 4 cores
library(BiocParallel)
register(SnowParam(8))

# load data
se <- read_rds('output/gene_expression.rds')

dat <- read.csv('output/gene_expression_averages.csv')
mat <- as.matrix(t(dat[, -1]))
rownames(mat) <- colnames(dat)[-1]
colnames(mat) <- dat[, 1]

## the tissue variable
g <- se$tissue[match(colnames(mat), se$cell_id)]

## the archetype variable
distances <- read_csv('output/cell_distances.csv')
gd <- spread(distances, archetype, dss) %>%
  column_to_rownames('cell_id') %>%
  apply(1, which.min)

all(names(gd) == colnames(mat))

df <- data.frame(tissue = g, archetype = as.factor(gd))

form <- ~ (1|tissue) + (1|archetype)

fit <- fitExtractVarPartModel(mat, form, df)
vp <- sortCols(fit)

write_csv(vp, 'output/variance_partition.csv')
