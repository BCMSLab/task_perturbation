# load libraries
library(tidyverse)

# load data
se <- read_rds('output/gene_expression.rds')

dat <- read.csv('output/gene_expression_averages.csv')
mat <- as.matrix(t(dat[, -1]))
rownames(mat) <- colnames(dat)[-1]
colnames(mat) <- dat[, 1]


mds <- cmdscale(dist(t(mat)))

## the tissue variable
g <- se$tissue[match(colnames(mat), se$cell_id)]

## the archetype variable
distances <- read_csv('output/cell_distances.csv')
gd <- spread(distances, archetype, dss) %>%
  column_to_rownames('cell_id') %>%
  apply(1, which.min)

all(names(gd) == colnames(mat))

mds <- as.data.frame(mds) %>%
  setNames(c('D1', 'D2')) %>%
  rownames_to_column('cell_id') %>%
  mutate(tissue = g, archetype = gd)

write_csv(mds, 'output/dimension_scaling.csv')

