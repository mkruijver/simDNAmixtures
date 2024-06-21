# hard coded constants related to the AMEL locus
AMEL_XX_unpacked <- c("X", "X")
AMEL_XX_packed <- "X,X"

AMEL_XY_unpacked <- c("X", "Y")
AMEL_XY_packed <- "X,Y"

dummy_pedigree <- list(ID = numeric(), FIDX = numeric(),
                       MIDX = numeric(), SEX = numeric())
class(dummy_pedigree) <- "ped"
