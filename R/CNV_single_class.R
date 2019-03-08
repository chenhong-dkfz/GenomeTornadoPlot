#' single CNV
#'
#' @export


CNV_single = setClass("CNV_single",
                      slots = list(
                        name = "character",
                        matrix = "GenomicRanges",
                        gene_name = "character"
                      ))
