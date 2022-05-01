fn_kit_data <- "data-raw/kit_data.txt"

x <- readr::read_tsv(fn_kit_data, col_types = "cccddddcddddccll")
kits <- c(list(all_kits = x), split(x, x$Panel))

usethis::use_data(kits, overwrite = TRUE, compress = 'xz')
