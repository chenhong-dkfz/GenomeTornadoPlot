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
#' @param type character, "dup" for all amplifications, "del" for all deletions. default none.
#' @examples
#' sdt <- MakeData(CNV=chr17,gene_name_1 = input_gene_1,score.type = "del")
#' where sdt is the intermediate file for TornadoPlot function
#'
#' @export


MakeData <- function(CNV,
                     gene_name_1,
                     gene_name_2,
                     score.type,max.length,score.method,cohort_thredshold,
                     gene_score_1,gene_score_2){

  data("genes",package = "GenomeTornadoPlot")

  if(missing(max.length)){
    max.length = 10000000 # 1e7
  }

  if(missing(cohort_thredshold)){
    cohort_thredshold <- 0.1
  }

  if(missing(gene_score_1)){
    gene_score_1 <- 9999
  }

  if(missing(gene_score_2)){
    gene_score_2 <- 9999
  }

  self_score <- TRUE
  if(gene_score_1==9999){self_score=FALSE}

  gene_coordinates = genes


  if(missing(gene_name_2)){   # in case there is only one gene of interests

    idx_gene_1 <- which(gene_coordinates$gene==gene_name_1)[1]
    if(gene_coordinates[idx_gene_1,"chromosome"]=="X"|gene_coordinates[idx_gene_1,"chromosome"]=="Y"|gene_coordinates[idx_gene_1,"chromosome"]=="M"){
      chrom <- as.character(gene_coordinates[idx_gene_1,"chromosome"])
    }else{
      chrom <- as.numeric(as.character(gene_coordinates[idx_gene_1,"chromosome"]))
    }
    start_1 <- as.numeric(as.character(gene_coordinates[idx_gene_1,"start"]))
    end_1 <- as.numeric(as.character(gene_coordinates[idx_gene_1,"end"]))


    if(missing(score.type)||score.type=="none"){
      score.type = "none"
    }else if(score.type=="dup"){
      score.type = "dup"
    }else{
      socre.type = "del"
    }  # so far score.type is only del#

    if(missing(score.method)){
      score.method = "normal"
    }

    CNV1 <- cbind(CNV,length=CNV$End-CNV$Start)
    CNV1 <- CNV1[CNV1$length >= 0,]
    CNV1 <- CNV1[CNV1$length <= max.length,]
    CNV1 <- CNV1[CNV1$Gene == gene_name_1,]
    CNV1 <- makeGRangesFromDataFrame(CNV1 , keep.extra.columns = TRUE)
    gene.position.1 <- GRanges(seqnames =Rle(chrom) , ranges=IRanges(start=start_1,end=end_1))
    CNV.gene1 <- subsetByOverlaps(CNV1,gene.position.1)


    if(score.type != "none"){
    fscore.cnv1 <- focallity.score.edge(gene_name=gene_name_1,cnv_file = CNV,filter = score.type,
                                        gene_coordinates = genes,method=score.method,
                                        max.length = max.length)
    }else{
      fscore.cnv1 <- 0
    }

    if(self_score==TRUE){fscore.cnv1<-gene_score_1}
    cnv_data <- new("CNV_single",name="CNV_test",matrix=CNV.gene1,gene_name=gene_name_1,gene_score=fscore.cnv1, t_gene_start = start_1,t_gene_end = end_1)


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


    CNV1 <- cbind(CNV,length=CNV$End-CNV$Start)
    CNV1 <- CNV1[CNV1$length >= 0,]
    CNV1 <- CNV1[CNV1$length <= max.length,]
    CNV1 <- CNV1[CNV1$Gene == gene_name_1,]

    CNV2 <- cbind(CNV,length=CNV$End-CNV$Start)
    CNV2 <- CNV2[CNV2$length >= 0,]
    CNV2 <- CNV2[CNV2$length <= max.length,]
    CNV2 <- CNV2[CNV2$Gene == gene_name_2,]

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


    if(missing(score.type)||score.type=="none"){
      score.type = "none"
    }else if(score.type=="dup"){
      score.type = "dup"
    }else{
      socre.type = "del"
    }  # so far score.type is only del#

    if(missing(score.method)){
      score.method = "normal"
    }

    if(score.type != "none"){
    fscore.cnv1 <- focallity.score.edge(gene_name_1,cnv_file = CNV,filter = score.type,
                                        gene_coordinates = genes,method=score.method,
                                        max.length = max.length)
    fscore.cnv2 <- focallity.score.edge(gene_name_2,cnv_file = CNV,filter = score.type,
                                        gene_coordinates = genes,method=score.method,
                                        max.length = max.length)
    }else{
      fscore.cnv1 <- 0
      fscore.cnv2 <- 0
    }

    if(self_score==TRUE){
      fscore.cnv1<-gene_score_1
      fscore.cnv2<-gene_score_2
    }

    cnv_data <- new("CNV_twin",name="Twin_Test",matrix_1=CNV1,
                    matrix_2=CNV2,gene_name_1=gene_name_1,gene_name_2=gene_name_2,
                    gene_score_1 = fscore.cnv1, gene_score_2 = fscore.cnv2,
                    t_gene_start_1 = start_1,t_gene_end_1 = end_1,
                    t_gene_start_2 = start_2,t_gene_end_2 = end_2,
                    max.length = max.length)
  }
  return(cnv_data)
}


