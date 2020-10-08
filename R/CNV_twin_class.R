CNV_twin = setClass("CNV_twin",
                    slots = list(
                      name = "character",
                      matrix_1 = "GenomicRanges",
                      matrix_2 = "GenomicRanges",
                      gene_name_1 = "character",
                      gene_name_2 = "character",
                      gene_score_1 = "numeric",
                      gene_score_2 = "numeric",
                      t_gene_start_1 = "numeric",
                      t_gene_end_1 = "numeric",
                      t_gene_start_2 = "numeric",
                      t_gene_end_2 = "numeric"
                    ))
