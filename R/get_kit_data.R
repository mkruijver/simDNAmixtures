#' Get kit data such as loci, dye, sizes for most standard kits
#'
#' @return Object of class list
#' @examples
#' # obtain the markers and dye colours for GlobalFiler
#' kit_data <- get_kit_data()
#' gf <- kit_data$GlobalFiler_Panel_v1
#' gf_markers <- unique(gf$Marker)
#' gf_colour_by_marker <- setNames(gf$Color[match(gf_markers, gf$Marker)], gf_markers)
#' @export
get_kit_data <- function() {
  fn_kit_data <- system.file("extdata/kit_data.txt", package = "SimMixDNA")

  x <- readr::read_tsv(fn_kit_data, col_types = "cccddddcddddccll")
  kit_data <- c(list(all_kits = x), split(x, x$Panel))

  kit_data
}
