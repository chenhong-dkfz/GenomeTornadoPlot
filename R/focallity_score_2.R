# input following data and parameters
# gene name, methods, filters, cnv file and gene coordinates

focallity.score.edge <- function(gene_name,method,filter,cnv_file,gene_coordinates,max.length){
  if(missing(method)){method="edge"}
    # if gene_1 is not unique in list, pick the first one
    # check gene_1 is the gene at the end of chromosome
    idx_gene_1 <- which(gene_coordinates$gene==gene_name)[1]
    gene_1_edge <- 0
    if(gene_coordinates[idx_gene_1,"chromosome"]=="X"|gene_coordinates[idx_gene_1,"chromosome"]=="Y"|gene_coordinates[idx_gene_1,"chromosome"]=="M"){
      gene_1_chr <- as.character(gene_coordinates[idx_gene_1,"chromosome"])
    }else{
      gene_1_chr <- as.numeric(as.character(gene_coordinates[idx_gene_1,"chromosome"]))
    }
    gene_1_start <- as.numeric(as.character(gene_coordinates[idx_gene_1,"start"]))
    gene_1_end <- as.numeric(as.character(gene_coordinates[idx_gene_1,"end"]))

    if(idx_gene_1==1){
      gene_1_edge <- 1
    }else if(idx_gene_1==nrow(gencode.v19.genes)){
      gene_1_edge <- 2
    }else if(as.numeric(gene_coordinates[idx_gene_1,"chromosome"]) != as.numeric(gene_coordinates[idx_gene_1-1,"chromosome"])){
      gene_1_edge <- 1
    }else if(as.numeric(gene_coordinates[idx_gene_1,"chromosome"]) != as.numeric(gene_coordinates[idx_gene_1+1,"chromosome"])){
      gene_1_edge <- 2
    }

    if(gene_1_edge==0){
      gene_1_left <-  gene_coordinates[idx_gene_1-1,"gene"]
      gene_1_left_start <- gene_coordinates[idx_gene_1-1,"start"]
      gene_1_left_end <-  gene_coordinates[idx_gene_1-1,"end"]
      gene_1_right <-  gene_coordinates[idx_gene_1-1,"gene"]
      gene_1_right_start <-  gene_coordinates[idx_gene_1+1,"start"]
      gene_1_right_end <-  gene_coordinates[idx_gene_1+1,"end"]
    }else if(gene_1_edge==1){
      gene_1_right <-  gene_coordinates[idx_gene_1+1,"gene"]
      gene_1_right_start <-  gene_coordinates[idx_gene_1+1,"start"]
      gene_1_right_end <-  gene_coordinates[idx_gene_1+1,"end"]
    }else if(gene_1_edge==2){
      gene_1_left <-  gene_coordinates[idx_gene_1-1,"gene"]
      gene_1_left_start <-  gene_coordinates[idx_gene_1-1,"start"]
      gene_1_left_end <-  gene_coordinates[idx_gene_1-1,"end"]
    }
    if(missing(max.length)){
      max.length = 10000000
    }

    if(missing(filter)){
      print("cnv scoring type missing, deletion as default")
      cnv_file <- cnv_file[cnv_file$Score<2,]
    }else if(type=="dup"){
      cnv_file <- cnv_file[cnv_file$Score>2,]
    }else if(type=="del"){
      cnv_file <- cnv_file[cnv_file$Score<2,]
    }

    CNV1 <- cbind(cnv_file,length=cnv_file$End-cnv_file$Start)
    CNV1 <- CNV1[CNV1$length >= 0,]
    CNV1 <- CNV1[CNV1$length <= as.numeric(as.character(max.length)),]
    CNV1 <- makeGRangesFromDataFrame(CNV1 , keep.extra.columns = TRUE)

    gene.position.1 <- GRanges(seqnames =Rle(gene_1_chr) , ranges=IRanges(start=gene_1_start,end=gene_1_end))
    CNV.gene1 <- subsetByOverlaps(CNV1,gene.position.1)
    if(gene_1_edge != 1){
      gene.position.l <- GRanges(seqnames =Rle(gene_1_chr) , ranges=IRanges(start=gene_1_left_start,end=gene_1_left_end))
      CNV.gene.l <- subsetByOverlaps(CNV1,gene.position.l)
      }
    if(gene_1_edge != 2){
      gene.position.r <- GRanges(seqnames =Rle(gene_1_chr) , ranges=IRanges(start=gene_1_right_start,end=gene_1_right_end))
      CNV.gene.r <- subsetByOverlaps(CNV1,gene.position.r)
    }

    starts_1 <- data.frame(CNV.gene1@ranges)$start
    ends_1 <- data.frame(CNV.gene1@ranges)$end
    scores_1 <- data.frame(CNV.gene1$Score)
    #max.length1 <- max(abs(ends_1-starts_1))
    max.length1 <- max.length
    focal.scores1 <- (log10(max.length1)-log10(ends_1-starts_1))*(scores_1+1)
    f.score1 <- sum(focal.scores1)
    if(!exists("CNV.gene.r")){
      starts_l <- data.frame(CNV.gene.l@ranges)$start
      ends_l <- data.frame(CNV.gene.l@ranges)$end
      scores_l <- data.frame(CNV.gene.l$Score)
      max.length_l <- max.length
      # max.length_l <- max(abs(ends_l-starts_l))
      focal.scores_l <- (log10(max.length_l)-log10(ends_l-starts_l))*(scores_l+1)
      f.score_d <- sum(focal.scores_l)
    }else if(!exists("CNV.gene.l")){
      starts_r <- data.frame(CNV.gene.r@ranges)$start
      ends_r <- data.frame(CNV.gene.r@ranges)$end
      scores_r <- data.frame(CNV.gene.r$Score)
      #max.length_r <- max(abs(ends_r-starts_r))
      max.length_r <- max.length
      focal.scores_r <- (log10(max.length_r)-log10(ends_r-starts_r))*(scores_r+1)
      f.score_d <- sum(focal.scores_r)
    }else{
      starts_l <- data.frame(CNV.gene.l@ranges)$start
      ends_l <- data.frame(CNV.gene.l@ranges)$end
      scores_l <- data.frame(CNV.gene.l$Score)
      #max.length_l <- max(abs(ends_l-starts_l))
      max.length_l <- max.length
      focal.scores_l <- (log10(max.length_l)-log10(ends_l-starts_l))*(scores_l+1)
      starts_r <- data.frame(CNV.gene.r@ranges)$start
      ends_r <- data.frame(CNV.gene.r@ranges)$end
      scores_r <- data.frame(CNV.gene.r$Score)
      #max.length_r <- max(abs(ends_r-starts_r))
      max.length_r <- max.length
      focal.scores_r <- (log10(max.length_r)-log10(ends_r-starts_r))*(scores_r+1)
      f.score_d <- 0.5*(sum(focal.scores_l)+sum(focal.scores_r))
    }
    if(method=="edge")
    {
      f.score <- f.score1 - f.score_d
    }else{
      f.score <- f.score1
    }

  f.score <- round(f.score,2)
  return(f.score)
}
