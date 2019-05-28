#' Generate input data for tornado plot
#'
#' @param CNV_1 data frame with six columns,Chromosome,Start,End,Score,Gene,Cohort,PID
#' @param CNV_2 data frame with six columns,Chromosome,Start,End,Score,Gene,Cohort,PID
#' @param gene_name_1 character, the name of first gene
#' @param gene_name_2 character, the name of second gene
#' @param start_1 numeric, start point of the first gene
#' @param start_2 numeric, start point of the second gene
#' @param end_1 numeric, end point of the first gene
#' @param end_2 numeric, end point of the second gene
#' @param chrom character, chromosome id of both genes
#' @param type character, "dup" for all duplications, "del" for all deletions. default none.
#'
#' @export


MakeData <- function(CNV_1,
                     gene_name_1,start_1,end_1,
                     CNV_2,
                     gene_name_2,start_2,end_2,
                     chrom,type){

  if(missing(CNV_2)){
    if(missing(type)){
      print("filter disabled!")
    }else if(type=="dup"){
      CNV_1 <- CNV_1[CNV_1$Score>2,]
    }else if(type=="del"){
      CNV_1 <- CNV_1[CNV_1$Score<2,]
    }
    CNV1 <- cbind(CNV_1,length=CNV_1$End-CNV_1$Start)
    CNV1 <- CNV1[CNV1$length >= 0,]
    CNV1 <- makeGRangesFromDataFrame(CNV1 , keep.extra.columns = TRUE)
    gene.position.1 <- GRanges(seqnames =Rle(chrom) , ranges=IRanges(start=start_1,end=end_1))
    CNV.gene1 <- subsetByOverlaps(CNV1,gene.position.1)
    cnv_data <- new("CNV_single",name="CNV_test",matrix=CNV.gene1,gene_name=gene_name_1)

  }else{
    if(missing(type)){
      print("filter disabled!")
    }else if(type=="dup"){
      CNV_1 <- CNV_1[CNV_1$Score>2,]
      CNV_2 <- CNV_2[CNV_2$Score>2,]
    }else if(type=="del"){
      CNV_1 <- CNV_1[CNV_1$Score<2,]
      CNV_2 <- CNV_2[CNV_2$Score<2,]
    }

    CNV1 <- cbind(CNV_1,length=CNV_1$End-CNV_1$Start)
    CNV1 <- CNV1[CNV1$length >= 0,]
    CNV1 <- makeGRangesFromDataFrame(CNV1 , keep.extra.columns = TRUE)
    gene.position <- GRanges(seqnames =Rle(chrom) , ranges=IRanges(start=start_1,end=end_1))
    CNV.gene1 <- subsetByOverlaps(CNV1,gene.position)

    CNV2 <- cbind(CNV_2,length=CNV_2$End-CNV_2$Start)
    CNV2 <- CNV2[CNV2$length >= 0,]
    CNV2 <- makeGRangesFromDataFrame(CNV2 , keep.extra.columns = TRUE)
    gene.position.2 <- GRanges(seqnames =Rle(chrom) , ranges=IRanges(start=start_2,end=end_2))
    CNV.gene2 <- subsetByOverlaps(CNV2,gene.position.2)
    cnv_data <- new("CNV_twin",name="Twin_Test",matrix_1=CNV.gene1,
                    matrix_2=CNV.gene2,gene_name_1=gene_name_1,gene_name_2=gene_name_2)
  }
  return(cnv_data)
}
