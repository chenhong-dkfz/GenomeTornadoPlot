#' Generate input data for tornado plot
#'
#' @param CNV data frame with six columns,Chromosome,Start,End,Score,Gene,Cohort,PID
#' @param CNV data frame with six columns,Chromosome,Start,End,Score,Gene,Cohort,PID
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


MakeData <- function(CNV,
                     gene_name_1,
                     gene_name_2,
                     score.type,max.length,score.method,cohort_thredshold){

  data("genes",package = "tornado.test.1")

  if(missing(max.length)){
    max.length = 10000000
  }

  if(missing(cohort_thredshold)){
    cohort_thredshold <- 0.1
  }


  gene_coordinates = genes

  print(nrow(gene_coordinates))


  if(missing(gene_name_2)){   # in case there is only one gene of interests

    idx_gene_1 <- which(gene_coordinates$gene==gene_name_1)[1]
    if(gene_coordinates[idx_gene_1,"chromosome"]=="X"|gene_coordinates[idx_gene_1,"chromosome"]=="Y"|gene_coordinates[idx_gene_1,"chromosome"]=="M"){
      chrom <- as.character(gene_coordinates[idx_gene_1,"chromosome"])
    }else{
      chrom <- as.numeric(as.character(gene_coordinates[idx_gene_1,"chromosome"]))
    }
    start_1 <- as.numeric(as.character(gene_coordinates[idx_gene_1,"start"]))
    end_1 <- as.numeric(as.character(gene_coordinates[idx_gene_1,"end"]))


    # if(missing(type)){
    #   print("filter disabled!")
    # }else if(type=="dup"){
    #   CNV <- CNV[CNV$Score>2,]
    # }else if(type=="del"){
    #   CNV <- CNV[CNV$Score<2,]
    # }

    if(missing(score.type)){
      score.type = "del"
    }else if(score.type=="dup"){
      score.type = "dup"
    }else{
      socre.type = "del"
    }

    if(missing(score.method)){
      score.method = "edge"
    }

    CNV1 <- cbind(CNV,length=CNV$End-CNV$Start)
    CNV1 <- CNV1[CNV1$length >= 0,]
    CNV1 <- CNV1[CNV1$length <= max.length,]
    CNV1 <- makeGRangesFromDataFrame(CNV1 , keep.extra.columns = TRUE)
    gene.position.1 <- GRanges(seqnames =Rle(chrom) , ranges=IRanges(start=start_1,end=end_1))
    CNV.gene1 <- subsetByOverlaps(CNV1,gene.position.1)

    fscore.cnv1 <- focallity.score.edge(gene_name_1,cnv_file = CNV,filter = score.type,
                                        gene_coordinates = gencode.v19.genes,method=score.method)
    cnv_data <- new("CNV_single",name="CNV_test",matrix=CNV.gene1,gene_name=gene_name_1,gene_score=fscore.cnv1)


  }else{   # in case there are two genes of interests

    idx_gene_1 <- which(gene_coordinates$gene==gene_name_1)[1]
    if(gene_coordinates[idx_gene_1,"chromosome"]=="X"|gene_coordinates[idx_gene_1,"chromosome"]=="Y"|gene_coordinates[idx_gene_1,"chromosome"]=="M"){
      chrom <- as.character(gene_coordinates[idx_gene_1,"chromosome"])
    }else{
      chrom <- as.numeric(as.character(gene_coordinates[idx_gene_1,"chromosome"]))
    }
    start_1 <- as.numeric(as.character(gene_coordinates[idx_gene_1,"start"]))
    end_1 <- as.numeric(as.character(gene_coordinates[idx_gene_1,"end"]))

    idx_gene_2 <- which(gene_coordinates$gene==gene_name_2)[1]
    if(gene_coordinates[idx_gene_2,"chromosome"]=="X"|gene_coordinates[idx_gene_2,"chromosome"]=="Y"|gene_coordinates[idx_gene_2,"chromosome"]=="M"){
      chrom <- as.character(gene_coordinates[idx_gene_2,"chromosome"])
    }else{
      chrom <- as.numeric(as.character(gene_coordinates[idx_gene_2,"chromosome"]))
    }
    start_2 <- as.numeric(as.character(gene_coordinates[idx_gene_2,"start"]))
    end_2 <- as.numeric(as.character(gene_coordinates[idx_gene_2,"end"]))


    if(missing(type)){
      print("filter disabled!")
    }else if(type=="dup"){
      CNV <- CNV[CNV$Score>2,]
    }else if(type=="del"){
      CNV <- CNV[CNV$Score<2,]
    }

    CNV1 <- cbind(CNV,length=CNV$End-CNV$Start)
    CNV1 <- CNV1[CNV1$length >= 0,]
    CNV1 <- CNV1[CNV1$length <= max.length,]

    CNV2 <- cbind(CNV,length=CNV$End-CNV$Start)
    CNV2 <- CNV2[CNV2$length >= 0,]
    CNV2 <- CNV2[CNV2$length <= max.length,]

    CNV1 <- makeGRangesFromDataFrame(CNV1 , keep.extra.columns = TRUE)
    gene.position <- GRanges(seqnames =Rle(chrom) , ranges=IRanges(start=start_1,end=end_1))
    CNV.gene1 <- subsetByOverlaps(CNV1,gene.position)
    CNV1 <- CNV.gene1


    CNV2 <- makeGRangesFromDataFrame(CNV2 , keep.extra.columns = TRUE)
    gene.position.2 <- GRanges(seqnames =Rle(chrom) , ranges=IRanges(start=start_2,end=end_2))
    CNV.gene2 <- subsetByOverlaps(CNV2,gene.position.2)
    CNV2 <- CNV.gene2

    CNV1$comp <- paste(CNV1@ranges,CNV1$PID,sep="_")
    CNV2$comp <- paste(CNV2@ranges,CNV2$PID,sep="_")
    CNV1$rep <- CNV1$comp %in% CNV2$comp
    CNV2$rep <- CNV2$comp %in% CNV1$comp

    if(nrow(data.frame(CNV1[CNV1$rep=="TRUE",]))!=0){CNV1[CNV1$rep=="TRUE",]$rep<- "C1"}

    if(nrow(data.frame(CNV2[CNV2$rep=="TRUE",]))!=0){CNV2[CNV2$rep=="TRUE",]$rep<- "C2"}

    if(nrow(data.frame(CNV1[CNV1$rep=="FALSE",]))!=0){CNV1[CNV1$rep=="FALSE",]$rep<- "U1"}

    if(nrow(data.frame(CNV2[CNV2$rep=="FALSE",]))!=0){CNV2[CNV2$rep=="FALSE",]$rep<- "U2"}


    fscore.cnv1 <- focallity.score.edge(gene_name_1,cnv_file = CNV,
                                        gene_coordinates = genes,method=score.method)
    fscore.cnv2 <- focallity.score.edge(gene_name_2,cnv_file = CNV,
                                        gene_coordinates = genes,method=score.method)
    cnv_data <- new("CNV_twin",name="Twin_Test",matrix_1=CNV1,
                    matrix_2=CNV2,gene_name_1=gene_name_1,gene_name_2=gene_name_2,
                    gene_score_1 = fscore.cnv1, gene_score_2 = fscore.cnv2)
  }
  return(cnv_data)
}


