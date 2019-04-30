#' @title splitBin
#' @description split transcriptome to sliding windows with step
#' @param gtf The gtf format gene annotation file
#' @param binSize Txdb file. Could be NULL when gtf is avaliable
#' @param step Step of Sliding window
#' @param threads Threads used for parallel running
#' @import data.table
#' @export


#
sub.win.gtf <- function(gtf,
                        binSize = 100,
                        step = 10,
                        threads = 4
                        )
{

  if(is.character(gtf)) {
    txdb <- makeTxDbFromGFF(gtf,format="gtf")
  } else {
    txdb <- gtf
  }
  # drop antisenses genes
  keepGene <- filterGene(txdb)
  exons = reduce(exonsBy(txdb,by="gene"))
  anno <- exons[names(exons) %in% names(keepGene),]

  registerDoParallel( cores = threads)
  newGTF <- foreach(ge = names(anno)) %dopar%{

    df.gene = as.data.table(anno[ge][[1]])
    exon.no = nrow(df.gene)
    exon.len = sum(df.gene$width)
    ## drop exons smaller than 200bp
    if(sum(df.gene$width)  < 200) {return(NULL)}

    if(exon.len <= binSize){
      slidingStart = 1
      slidingEnd = exon.len
    } else {
      slidingStart =round(seq(from = 1, to = exon.len-binSize, by = step))
      slidingEnd = c((slidingStart+binSize)[-length(slidingStart)],exon.len)
    }

    site = 1:exon.len; start = 1
    for (no.E in 1:exon.no) {
      site[start:(start+df.gene$width[no.E]-1)] <- c(df.gene$start[no.E]:df.gene$end[no.E])
      start = start+df.gene$width[no.E]
    }
    x1 <- data.table(anno[ge][[1]]@seqnames@values,site[slidingStart],site[slidingEnd],bin=1:length(slidingStart))
    colnames(x1) <- c("seqnames","start","end","bin")

    res <- vector("list",nrow(x1))
    for (n in 1:length(res)) {
      res[[n]] <- n
    }
    aa <- lapply(res,covertNoIntron,x1=x1,df.gene=df.gene)
    aa <- rbindlist(aa); aa <- aa[order(aa$bin)]
    data.table(aa[,1],".","bin",aa[,2:3],".",anno[ge][[1]]@strand@values,".",as.data.table(paste('bin_id "',paste(geneName,aa[,4][[1]],sep = "-"),'";',sep="")))

  }

  return(rbindlist(newGTF))

}

covertNoIntron <- function(line,x1=x1,df.gene=df.gene) {

  num1 <- which((x1[line,2][[1]]<=df.gene$end)&(x1[line,2][[1]]>=df.gene$start))
  num2 <- which((x1[line,3][[1]]<=df.gene$end)&(x1[line,3][[1]]>=df.gene$start))

  if(num2==num1) {
    return(x1[line,])
  } else  {
    a <- x1[line,][rep(1,each=(num2-num1+1)),]
    a[1,3] <- df.gene[num1[1],3][[1]]
    a[2:(num2-num1+1),2] <- df.gene[c(num1:num2)[2:(num2-num1+1)],2][[1]]

    if(num2-num1 > 1) {
      a[2:(num2-num1),3] <- df.gene[c(num1:num2)[2:(num2-num1)],3]
    }
  }
  return(a)
}
