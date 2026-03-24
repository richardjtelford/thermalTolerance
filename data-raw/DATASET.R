## code to prepare `htdata` dataset goes here

htdata <- read.csv("data-raw/Sample_FvFm_data.csv")
colnames(htdata) <- tolower(colnames(htdata))

usethis::use_data(htdata, overwrite = TRUE)
