# load libraries ----
library(SummarizedExperiment)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(archetypes)
library(circlize)
library(cowplot)
library(variancePartition)

# ggplot theme ----
my_theme <- theme(panel.grid = element_blank(),
                  panel.border = element_rect(fill = NA, color = 'black', size = 1),
                  panel.background = element_blank(),
                  axis.title.x.top = element_blank(),
                  axis.text.x.top = element_blank(),
                  axis.title.y.right = element_blank(),
                  axis.text.y.right = element_blank(),
                  legend.background = element_blank(),
                  axis.ticks.length=unit(2, "mm"))

# load data ----
cancer_cells_genes_average <- read_rds('data_clean/cancer_cells_genes_average.rds')
cancer_cells_scores <- read_rds('data_clean/cancer_cells_scores.rds')
breast_cells_growth <- read_rds('data_clean/breast_cells_growth.rds')
cell_tissue <- read_rds('data_clean/cell_tissue.rds')

cell_archetypes <- read_rds('data_clean/cell_archetypes.rds')
archetypes_6 <- read_rds('data_clean/archetypes_6.rds')
arc_cell_distance <- read_rds('data_clean/arc_cell_distances.rds')
enrich_treatment <- read_rds('data_clean/enrich_treatment.rds')

# correlation average expression with fold-change ----
## extract the names of the cell lines
cells <- str_split(names(cancer_cells_scores), '_', simplify = TRUE)[, 1]
ind <- unique(cells)
names(ind) <- ind

## loop over the cells and extract fold-change
mats <- map(ind, function(c) {
  # subset the the samples for each cell
  l <- cancer_cells_scores[cells == c]
  
  # join the data for each cell
  df <- Reduce(left_join, l)
  
  # remove control samples and other variables
  # remove NA and duplicates
  # transform to a matrix
  df %>%
    dplyr::select(-gene, -cell_id, -DMSO) %>%
    filter(!duplicated(symbol)) %>%
    filter(!is.na(symbol)) %>%
    column_to_rownames('symbol') %>%
    as.matrix()
})

# remove cells with only one drug
mats <- mats[map(mats, ncol) > 1]

# extract common cell lines
ind <- intersect(names(mats), rownames(cancer_cells_genes_average))
names(ind) <- ind

# for each cell calculate the correaltions
corrs <- map_df(ind, function(c) {
  # fold change
  x <- mats[[c]]
  
  # average expression
  y <- unlist(cancer_cells_genes_average[c,])
  
  # calculate the correlation between variables
  # transform to tibble
  sym <- intersect(rownames(x), names(y))
  as_tibble(melt(cor(x[sym,], y[sym], method = 'spearman')))
}, .id = 'cell_id')

## to use in the report
length(unique(corrs$cell_id))
length(unique(corrs$Var1))

set.seed(123)
ks.test(corrs$value,
        'rnorm',
        alternative = 'less')

drug_corrs <- corrs %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = NA, color = 'black', bins = 20) +
  geom_vline(xintercept = 0, lty = 2, color = 'red') +
  labs(x = "Spearman's rank correlation",
       y = 'Number of drugs') +
  scale_x_continuous(sec.axis = dup_axis(), limits = c(-.32, .32)) +
  scale_y_continuous(sec.axis = dup_axis()) +
  my_theme

group_corrs <- list('Drug' = group_by(corrs, Var1) %>% summarise(value = mean(value)) %>% pull(value),
     'Cell Line' = group_by(corrs, cell_id) %>% summarise(value = mean(value)) %>% pull(value)) %>%
  melt()

set.seed(123)
split(group_corrs$value, group_corrs$L1) %>%
  map(ks.test, y = 'rnorm', alternative = 'less')

average_corrs <- group_corrs %>%
  ggplot(aes(x = value, fill = L1)) +
  geom_histogram(bins = 40, alpha = .5) +
  geom_vline(xintercept = 0, lty = 2, color = 'red') +
  scale_y_continuous(sec.axis = dup_axis(),
                     breaks = seq(0, 10, 2),
                     limits = c(0, 7)) +
  scale_x_continuous(sec.axis = dup_axis(),
                     limits = c(-.25, .25),
                     breaks = c(-0.2, 0, 0.2)) +
  labs(x = "Spearman's rank correlation",
       y = 'Number of drugs/cell lines',
       fill = '') +
  my_theme +
  theme(legend.position = c(.75, .85))

plot_grid(drug_corrs,
          average_corrs,
          scale = .9, 
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/expression_fc_correlation.png',
         height = 2.8, width = 5.6)

# multidimensional scaling ----
## calculate mds
mds <- cmdscale(dist(cancer_cells_genes_average))

## the tissue variable
g <- cell_tissue$tissue[match(rownames(mds), cell_tissue$cell_id)]

## the archetype variable
gd <- apply(arc_cell_distance, 1, which.min)

df <- data.frame(tissue = g, archetype = as.factor(gd))
head(df)
form <- ~ (1| tissue) + (1|archetype)

fit <- fitExtractVarPartModel(as.matrix(t(cancer_cells_genes_average)),
                              form,
                              df)
vp <- sortCols(fit)

median(vp$tissue)
median(vp$archetype)

mds_tissue <- as_tibble(mds) %>%
  ggplot(aes(x = V1, y = V2, color = g)) +
  geom_point(shape = 'square') +
  labs(x = 'D1', y = 'D2') +
  scale_x_continuous(sec.axis = dup_axis(),
                     breaks = c(-60, -30, 0, 30, 60)) +
  scale_y_continuous(sec.axis = dup_axis(),
                     breaks = c(-60, -30, 0, 30, 60)) +
  my_theme +
  theme(legend.position = 'none')

mds_archetype <- as_tibble(mds) %>%
  ggplot(aes(x = V1, y = V2, color = as.factor(gd))) +
  geom_point(shape = 'square') +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = 'D1', y = 'D2') +
  scale_x_continuous(sec.axis = dup_axis(),
                     breaks = c(-60, -30, 0, 30, 60)) +
  scale_y_continuous(sec.axis = dup_axis(),
                     breaks = c(-60, -30, 0, 30, 60)) +
  my_theme +
  theme(legend.position = 'none')

length(unique(rownames(mds)))

plot_grid(mds_tissue,
          mds_archetype,
          scale = .9, 
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/mds.png',
         height = 2.8, width = 5.6)

# plot rss and archetypes ----
rss_df <- tibble(x = 1:10, y = rss(cell_archetypes)[, 1])

rss_plot <-  rss_df %>%
  ggplot(aes(x, y)) +
  geom_line() +
  geom_hline(yintercept = rss_df$y[rss_df$x == 6], lty = 2, color = 'red') +
  # scale_x_continuous(breaks = seq(2, 10, 2)) +
  # scale_y_continuous(breaks = seq(2, 8, 2), limits = c(2, 8)) +
  labs(x = 'Number of Archetypes',
       y = 'Residual Sum of Square Errors') +
  scale_x_continuous(sec.axis = dup_axis(), limits = c(1, 10), breaks = seq(2, 10, 2)) +
  scale_y_continuous(sec.axis = dup_axis(), limits = c(1, 10), breaks = seq(2, 10, 2)) +
  my_theme

## simplex
simp <- simplexplot(archetypes_6)
proj_h <- as_tibble(simp$proj_h)
proj_z <- as_tibble(simp$proj_z)

simp_plot <- ggplot() +
  geom_point(data = proj_h, aes(x = x, y = y), shape = 'square') +
  geom_text(data = proj_z, aes(x = x, y = y,
                                label = paste0('A', 1:6)),
             color = 'red') +
  scale_x_continuous(sec.axis = dup_axis(), limits = c(-8.5, 8.5)) +
  scale_y_continuous(sec.axis = dup_axis(), limits = c(-8.5, 8.5)) +
  my_theme

plot_grid(rss_plot,
          simp_plot,
          scale = .9, 
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/archeytpes_rss.png',
         height = 2.8, width = 5.6)

# distances ----
d <- arc_cell_distance
d[d > 120] <- max(d)
col_fun <- colorRamp2(c(max(d), min(d)),
                      c('white', 'black'))
colnames(d) <- paste0('A', colnames(d))

g <- cell_tissue$tissue[match(rownames(d), cell_tissue$cell_id)]

# Heatmap(d,
#         col = col_fun,
#         column_names_rot = 0,
#         show_heatmap_legend = FALSE,
#         show_row_dend = FALSE,
#         show_column_dend = FALSE,
#         row_split = g,
#         column_names_centered = TRUE)

d2 <- d
d2[d > 120] <- 0
d2[d < 120] <- 1

tibble(tissue = g,
       as.data.frame(d2) %>% rownames_to_column('cell_id')) %>%
  arrange(tissue, cell_id, A1, A2, A3, A4, A5, A6) %>%
  xtable::xtable()

# hallmark enrichment ----
# png(filename = 'output/manuscript/figures/hallmark_enrichment.png',
#     width = 5, height = 3, units = 'in', res = 300)
# enrich_archetypes %>%
#   filter(p.adjust < .2) %>%
#   mutate(archetype = paste0('A', archetype)) %>%
#   acast(ID ~ archetype, value.var = 'NES') %>%
#   Heatmap( column_names_rot = 0,
#            show_heatmap_legend = FALSE,
#            show_row_dend = FALSE,
#            show_column_dend = FALSE,
#            column_names_centered = TRUE)
# dev.off()

enrich_archetypes <- read_rds('data_clean/enrich_archetype.rds')

enrich_archetypes_filtered <- enrich_archetypes %>%
  filter(p.adjust < .2) %>%
  mutate(archetype = paste0('A', archetype))

enrich_archetypes_plot <- enrich_archetypes_filtered %>%
  ggplot(aes(x = archetype,
             y = ID,
             size = abs(NES),
             color = NES > 0)) +
  geom_point() +
  scale_size_continuous(breaks = c(1, 1.5, 2), limits = c(1, 2.2)) +
  labs(x = 'Archetype', y = '', size = 'NES') +
  my_theme

ggsave(plot = enrich_archetypes_plot,
         filename = 'output/manuscript/figures/hallmark_enrichment.png',
         width = 5.6, height = 2.8)

# correlation with drug enrichment ----

## merge cell distances and enrichment
enrich <- arc_cell_distance %>%
  melt() %>%
  as_tibble() %>%
  setNames(c('cell_id', 'archetype', 'distance')) %>%
  filter(distance < 120) %>%
  mutate(archetype = as.character(archetype)) %>%
  left_join( dplyr::select(enrich_archetypes, archetype, ID, archetype_es = NES)) %>%
  left_join(dplyr::select(enrich_treatment, cell_id, drug, ID, drug_es = NES)) %>%
  group_by(cell_id, archetype, distance, ID) %>%
  # top_n(5, abs(drug_es)) %>%
  summarise_at(vars(ends_with('es')), mean)

## make groups based on enrichment in archetypes
enrich_group <- enrich %>%
  ungroup() %>%
  mutate(group = cut(round(abs(drug_es), 1), 4,
                     labels = c('< .47', '(.47,.95]', '(.95,1.4]', '(1.4,1.9]'))) %>%
  group_by(group) %>%
  summarise(corr = cor(archetype_es, drug_es, use = 'complete'))

enrich_group_plot <- enrich_group %>%
  ggplot(aes(x = as.numeric(group), y = corr)) +
  geom_col(fill = NA, color = 'black', width = .7) +
  labs(x = 'Enrichment in archetype (absolute)',
       y = 'PCC; with drug response') +
  scale_x_continuous(breaks = 1:4,
                     labels = enrich_group$group,
                     sec.axis = dup_axis()) +
  scale_y_continuous(sec.axis = dup_axis()) +
  my_theme

# plot correlation between expression change and growth rate
correlation_group <- breast_cells_growth %>%
  select(cell_id = Cell.Name,
         Var1 = Small.Molecule.Name,
         conc = Perturbagen.Concentration,
         grmax = Normalized.Growth.Rate.Inhibition.Value) %>%
  mutate(Var1 = tolower(Var1)) %>%
  inner_join(corrs) %>%
  group_by(cell_id, Var1) %>%
  top_n(-1, grmax) %>%
  unique() %>%
  ungroup() %>%
  mutate(group = cut((value),
                     c(-.4, -.3,-.2, -.1, 0, .1, .2),
                     labels = c('< -.3', '(-.3,-.2]', '(-.2,-.1]', '(-.1,0]', '(0,.1]', '(.1,.2]')))

correlation_group_grmax <- correlation_group %>%
  ggplot(aes(x = as.numeric(group), y = grmax, group = group)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2, color = 'red') +
  labs(x = 'Expression Change Correlation',
       y = 'Max Growth Inhibition') +
  my_theme +
  theme(axis.text.x = element_text(size = 6.5)) +
  scale_x_continuous(breaks = 1:length(correlation_group$group),
                     labels = correlation_group$group,
                     sec.axis = dup_axis()) +
  scale_y_continuous(sec.axis = dup_axis()) +
  my_theme

plot_grid(enrich_group_plot,
          correlation_group_grmax,
          scale = .9, 
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/enrich_groups.png',
         height = 2.8, width = 5.6)
