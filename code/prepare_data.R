# load libraries
library(tidyverse)
library(reshape2)
library(slinky)
library(org.Hs.eg.db)

# # Breast Cancer Profiling Project, Gene Expression 1: Baseline mRNA sequencing on 35 breast cell lines - Dataset (ID:20348)
# url <- 'http://lincs.hms.harvard.edu/data/HMS_Dataset_20348.zip'
# download.file(url, destfile = 'data/HMS_Dataset_20348.zip')
# unzip('data/HMS_Dataset_20348.zip', exdir = 'data')
# 
# # Breast Cancer Profiling Project – Proteomics 1: 1 total proteome dataset for a 35-cell line breast cancer panel under basal conditions - Dataset (ID:20352)
# url <- 'http://lincs.hms.harvard.edu/data/HMS_Dataset_20352.zip'
# download.file(url, destfile = 'data/HMS_Dataset_20352.zip')
# unzip('data/HMS_Dataset_20352.zip', exdir = 'data')
# 
# # Growth rate-corrected (GR) dose-response metrics across a panel of 71 breast cancer cell lines treated with a library of small molecule and antibody perturbagens. Dataset 1 of 4: Relative cell counts and normalized growth rate inhibition values across technical replicates. - Dataset (ID:20268)
# url <- 'https://lincs.hms.harvard.edu/db/datasets/20268/results?search=&output_type=.csv'
# download.file(url, destfile = 'data/HMS_Dataset_20268.csv')

# gene expression
breast_cells_genes_df <- read_csv('data/HMS_Dataset_20348/HMS_Dataset_20348_DataFile.csv')
breast_cells_genes <- as.matrix(breast_cells_genes_df[, -1])
rownames(breast_cells_genes) <- breast_cells_genes_df$id
breast_cells_genes <- breast_cells_genes[!duplicated(rownames(breast_cells_genes)),]

write_rds(breast_cells_genes, 'data_clean/breast_cells_genes.rds')

rm(breast_cells_genes_df)

# proteoms
breast_cells_proteins_df <- read_csv('data/HMS_Dataset_20352/HMS_Dataset_20352_DataFile.csv')
breast_cells_proteins <- as.matrix(breast_cells_proteins_df[, -c(1,2)])
rownames(breast_cells_proteins) <- breast_cells_proteins_df$Gene_symbol
breast_cells_proteins <- breast_cells_proteins[!duplicated(rownames(breast_cells_proteins)),]

write_rds(breast_cells_proteins, 'data_clean/breast_cells_proteins.rds')

rm(breast_cells_proteins_df)

# growth rate
breast_cells_growth <- read_csv('data/HMS_Dataset_20268.csv', name_repair = 'universal')

write_rds(breast_cells_growth, 'data_clean/breast_cells_growth.rds')

# differential expression

# download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_inst_info.txt.gz',
#               destfile = 'data/GSE92742_Broad_LINCS_inst_info.txt.gz')
# aria2c -x 8 -s 8 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz

user_key <- httr::content(httr::GET("https://api.clue.io/temp_api_key"), as = "parsed")$user_key
gctx_file <- '~/workingon/LINPS/data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
gctx_info <- '~/workingon/LINPS/data/GSE92742_Broad_LINCS_inst_info.txt.gz'

sl <- Slinky(user_key, gctx_file, gctx_info)
md <- slinky::metadata(sl)

pert_iname <- unique(md$pert_iname[md$pert_iname %in% tolower(breast_cells_growth$Small.Molecule.Name)])
breast_cells <- unique(md$cell_id[md$cell_id %in% breast_cells_growth$Cell.Name])
cancer_cells <- unique(md$cell_id[md$pert_iname %in% pert_iname])

id <- which((md$pert_iname %in% c(pert_iname, 'DMSO')) & (md$cell_id %in% cancer_cells))
se <- as(sl[, id], 'SummarizedExperiment')

dir.create('scores')

# make sure this is correct!!

# for(cell in unique(se$cell_id)) {
#   lapply(unique(se$pert_iname),
#            function(drug) {
#              n <- sum(se$cell_id == cell & se$pert_iname == drug)
#              if (n > 5) {
#                ks_vecs <- diffexp(
#                  sl,
#                  treat = se[, se$cell_id == cell & se$pert_iname == drug],
#                  control = se[, se$cell_id == cell & se$pert_iname == 'DMSO'],
#                  method = 'ks',
#                  verbose = TRUE
#                )
#                
#                df <- melt(ks_vecs)
#                names(df) <- c('gene', drug)
#                df$cell_id <- cell
#                
#                write_csv(df, file.path('scores', paste0(cell, '_', drug, '.csv')))
#              }
#            })
# }

fls <- list.files('data_clean/scores/', full.names = TRUE)
names(fls) <- str_split(fls, '\\.|/', simplify = TRUE)[, 3]

id_symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                   rownames(se),
                                   'SYMBOL',
                                   'ENTREZID') %>%
  na.omit() %>%
  filter(!duplicated(SYMBOL)) %>%
  dplyr::select(gene = ENTREZID, symbol = SYMBOL) %>%
  as_tibble() %>%
  mutate(gene = as.integer(gene))

cancer_cells_scores <- map(fls, function(x) {
  df <- read_csv(x)
  left_join(df, id_symbol)
})

write_rds(cancer_cells_scores, 'data_clean/cancer_cells_scores.rds')

# cancer cells genes
id <- which((md$pert_iname == 'DMSO') & (md$cell_id %in% cancer_cells))
se <- as(sl[, id], 'SummarizedExperiment')

cell_tissue <- read_tsv('data/cellinfo_beta.txt')
cell_tissue <- cell_tissue[match(se$cell_id, cell_tissue$cell_iname),]
se$tissue <- cell_tissue$cell_lineage
se$subtype <- cell_tissue$subtype

id_symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                   rownames(se),
                                   'SYMBOL',
                                   'ENTREZID')

rownames(se) <- id_symbol$SYMBOL

write_rds(se, 'data_clean/cancer_cells_genes.rds')

# terms to genes
term_gene <- select(
  org.Hs.eg.db,
  rownames(se),
  'GO',
  'SYMBOL'
) %>%
  filter(ONTOLOGY == 'BP') %>%
  mutate(term = unlist(Term(GO.db::GOTERM[GO]), use.names = FALSE)) %>%
  dplyr::select(term, gene = SYMBOL) %>%
  unique() %>%
  as_tibble()

write_rds(term_gene, 'data_clean/term_gene.rds')

# hallmarks
# source: CHG: A Systematically Integrated Database of Cancer Hallmark Genes
# source: http://bio-bigdata.hrbmu.edu.cn/CHG/index.html
url <- 'http://bio-bigdata.hrbmu.edu.cn/CHG/download/Supplementary%20Table%202%20Genesets%20of%2010%20hallmarker.xlsx'

download.file(url, 'data/hallmarks_gene_sets.xlsx')
symbols <- keys(org.Hs.eg.db, 'SYMBOL')

hallmarks <- readxl::read_excel('data/hallmarks_gene_sets.xlsx',
                   skip = 2) %>%
  as.list() %>%
  map(function(x) {
    l <- x[!is.na(x)]
    l <- toupper(l)
    intersect(l, symbols)
  })

write_rds(hallmarks, 'data_clean/hallmarks_sets.rds')

