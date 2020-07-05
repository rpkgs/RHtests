## code to prepare `DATASET` dataset goes here
library(data.table)

PTtable <- fread("data-raw/PTtable.csv", skip = 2) %>% as.matrix()
PFtable <- fread("data-raw/PFtable.csv") %>% as.matrix()

# usethis::use_data(DATASET, overwrite = TRUE)

use_data(PFtable, overwrite = TRUE)
use_data(PTtable, overwrite = TRUE)
