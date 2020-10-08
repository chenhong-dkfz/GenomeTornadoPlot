#' single CNV
#'
#' @export


CNV_single = setClass("CNV_single",
                      slots = list(
                        name = "character",
                        matrix = "GenomicRanges",
                        gene_name = "character",
                        gene_score = "numeric",
                        t_gene_start = "numeric",
                        t_gene_end = "numeric"
                      ))
