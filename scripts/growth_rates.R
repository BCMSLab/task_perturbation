# load libraries
library(tidyverse)

# # Growth rate-corrected (GR) dose-response metrics across a panel of 71 breast cancer cell lines treated with a library of small molecule and antibody perturbagens. Dataset 1 of 4: Relative cell counts and normalized growth rate inhibition values across technical replicates. - Dataset (ID:20268)
# url <- 'https://lincs.hms.harvard.edu/db/datasets/20268/results?search=&output_type=.csv'
# download.file(url, destfile = 'data/HMS_Dataset_20268.csv')

# read
growth_rates <- read_csv('data/HMS_Dataset_20268.csv', name_repair = 'universal')

# tidy
growth_rates <- growth_rates %>%
  dplyr::select(cell_id = Cell.Name,
                drug = Small.Molecule.Name,
                dose = Perturbagen.Concentration,
                ndr = Nominal.Division.Rate,
                grmax = Normalized.Growth.Rate.Inhibition.Value) %>%
  mutate(drug = tolower(drug),
         cell_id = str_replace_all(cell_id, '-', ''))

write_csv(growth_rates, 'output/growth_rates.csv')
