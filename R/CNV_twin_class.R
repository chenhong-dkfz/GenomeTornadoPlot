CNV_twin = setClass("CNV_twin",
                    slots = list(
                      name = "character",
                      matrix_1 = "GenomicRanges",
                      matrix_2 = "GenomicRanges",
                      gene_name_1 = "character",
                      gene_name_2 = "character"
                    ))
