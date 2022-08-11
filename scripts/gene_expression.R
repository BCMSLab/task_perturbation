# Extract gene expression from LINCS

# load libraries
library(tidyverse)
library(slinky)
library(org.Hs.eg.db)

# download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_inst_info.txt.gz',
#               destfile = 'data/GSE92742_Broad_LINCS_inst_info.txt.gz')
# aria2c -x 8 -s 8 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz

# link to gctx files and key
user_key <- httr::content(httr::GET("https://api.clue.io/temp_api_key"), as = "parsed")$user_key
gctx_file <- 'data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
gctx_info <- 'data/GSE92742_Broad_LINCS_inst_info.txt.gz'

# load gctx file and metadata
sl <- Slinky(user_key, gctx_file, gctx_info)
md <- slinky::metadata(sl)

# load the names of drugs from growth_rates
growth_rates <- read_csv('output/growth_rates.csv')
breast_cell_drugs <- unique(growth_rates$drug)

# extract common ids 
pert_iname <- unique(md$pert_iname[md$pert_iname %in% tolower(breast_cell_drugs)])
cancer_cells <- unique(md$cell_id[md$pert_iname %in% pert_iname])


id <- which((md$pert_iname == 'DMSO') & (md$cell_id %in% cancer_cells))
se <- as(sl[, id], 'SummarizedExperiment')

# cell information
url <- "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt"
download.file(url, 'data/cellinfo_beta.txt')
cell_tissue <- read_tsv('data/cellinfo_beta.txt')

cell_tissue <- cell_tissue[match(se$cell_id, cell_tissue$cell_iname),]
se$tissue <- cell_tissue$cell_lineage
se$subtype <- cell_tissue$subtype

id_symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  rownames(se),
  'SYMBOL',
  'ENTREZID'
)

rownames(se) <- id_symbol$SYMBOL

write_rds(se, 'output/gene_expression.rds')
