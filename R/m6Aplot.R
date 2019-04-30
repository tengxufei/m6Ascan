#' @title m6Aplot
#' @description  plot peak distribution among elements of transcriptome and genome
#' @param species The species of your sample
#' @param gtf The gtf format gene annotation file, only used when species is empty
#' @param frac The length of single element
#' @param strand overlap gene with peak at same strand
#' @param peakFile Path to the peak file
#' @export

tryCatch(suppressPackageStartupMessages(require(data.table)),error=install.packages("data.table",repos = "http://cran.r-project.org"),finally=print("package 'data.table' is needed, it's downloading ...") )
tryCatch(suppressPackageStartupMessages(require(ggplot2)),error=install.packages("ggplot2",repos = "http://cran.r-project.org"),finally=print("package 'ggplot2' is needed, it's downloading ...") )

plot.peaks <- function(txdb=NULL,# only txdb like "TxDb.Mmusculus.UCSC.mm10.knownGene" is support
                       gtf, # gtf file used for gene extraction
                       frac = c(20,100,80),
                       strand = "all", # could be "all","positive" and "negative",
                       species = "mm10", # only used for mouse or human, "txdb" option is universe
                       peak)
{
  # get gene annotation file
  if(is.null(txdb)) {
    if(is.character(gtf)) {
      cat("Making txdb from gtf file ...","\n")
      txdb=makeTxDbFromGFF(gtf,format="gtf")
    } else {
      if(species == "mm10") {
        require(TxDb.Mmusculus.UCSC.mm10.knownGene)
        cat("Making txdb from UCSC mm10 annotation file ...","\n")
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
      } else if (species == "hg38") {
        require("TxDb.Hsapiens.UCSC.hg38.knownGene")
        cat("Making txdb from UCSC hg38 annotation file ...","\n")
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
      } else {
        cat("The txdb or gtf file is wrong, please check", "\n")
      }
    }
  }
  keepGene <- filterGene(txdb)
  peak <- as.data.frame(peak)
  colnames(peak) <- c("seqnames","start","end","name","strand")
  peak <- peak[peak$name %in% names(keepGene),]
  peak$seqnames <- factor(peak$seqnames)

  # check the chrmosome format
  if(any(grepl("chr",seqlevels(txdb)))) {
    if(!any(grepl("chr",levels(peak$seqnames))))  {
      levels(peak$seqnames) <- paste("chr",levels(peak$seqnames),sep="")
    }
  } else {
    if(any(grepl("chr",levels(peak$seqnames)))) {
      levels(peak$seqnames) <- levels(factor(gsub("chr","",peak$seqnames)))
    }
  }

  # make annotation file based on txdb
  cat("Making split elements ...","\n")
  annoM <- processTxdb(txdb=txdb,frac=frac)
  # calculate weights
  cat("Calculating the weight of each peak in each bin ...","\n")
  mat <- calPeak(annoM,peak,frac=frac)
  ylen <- max(mat$num)
  ggplot(mat,aes(x=win,y=num))+
    geom_area(col="skyblue",fill="skyblue",alpha=0.5,lwd=1.5)+
    ggtitle("m6A peak distribution on mRNA")+theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
    annotate("text", x = frac[1]/2, y = -(ylen/20), label = "5'UTR")+
    annotate("text", x = frac[1]+(frac[2]/2), y = -(ylen/20), label = "CDS")+
    annotate("text", x = sum(frac[1:2])+(frac[3]/2), y = -(ylen/20), label = "3'UTR")+
    geom_vline(xintercept=c(frac[1],sum(frac[1:2])), linetype="dotdash",col="grey",lwd=1) +
    annotate("rect", xmin = 1, xmax = frac[1], ymin = 0, ymax = -(ylen/35), alpha = .2, colour = "black") +
    annotate("rect", xmin = frac[1], xmax = sum(frac[1:2]), ymin = 0, ymax = -(ylen/35), alpha = .99, colour = "black") +
    annotate("rect", xmin = sum(frac[1:2]), xmax = sum(frac), ymin = 0, ymax = -(ylen/35), alpha = .2, colour = "black")
}

## drop antisenses genes
filterGene <- function(txdb)
{
  single.strand.genes <- genes(txdb, columns=c("tx_chrom", "tx_strand"))
  return(single.strand.genes)
}

## process Txdb to bin-split-format

processTxdb <- function(txdb,frac=c(20,100,80))
{
  trans <- transcripts(txdb)
  longestIsoform <- as.data.frame(reduce(trans)) #get the longest isoform
  longestIsoform <-longestIsoform[order(longestIsoform$seqnames,longestIsoform$start),]

  cds <- orderAndCom(cds(txdb),longestIsoform)
  utr5 <- orderAndCom(fiveUTRsByTranscript(txdb),longestIsoform)
  utr5 <- utr5[,c(1:5,8:9)]
  utr3 <- orderAndCom(threeUTRsByTranscript(txdb),longestIsoform)
  utr3 <- utr3[,c(1:5,8:9)]
  # when subjectHits is consistence, add the width together
  #utr5R <- tapply(utr5$width,utr5$subjectHits,function(x){round(frac[1]/sum(x),2)})
  #cdsR <- tapply(cds$width,cds$subjectHits,function(x){round(frac[2]/sum(x),2)})
  #utr3R <- tapply(utr3$width,utr3$subjectHits,function(x){round(frac[3]/sum(x),2)})

  utr5S <- tapply(utr5$width,utr5$subjectHits,sum)
  cdsS <- tapply(cds$width,cds$subjectHits,sum)
  utr3S <- tapply(utr3$width,utr3$subjectHits,sum)

  utr5M <- tapply(utr5$queryHits,utr5$subjectHits,min)
  cdsM <- tapply(cds$queryHits,cds$subjectHits,min)
  utr3M <- tapply(utr3$queryHits,utr3$subjectHits,min)

  result <- list('cds' = cds,'cdsM' = cdsM, 'cdsS' = cdsS,'utr5' = utr5,'utr5M' = utr5M, 'utr5S' = utr5S,'utr3' = utr3,'utr3M' = utr3M, 'utr3S' = utr3S)
  return(result)
}

##
calPeak <- function(annoM,peak,frac) {

  dataCDS <- overlapMatrix(annoM$cds,peak)
  dataUTR3 <- overlapMatrix(annoM$utr3,peak)
  dataUTR5 <- overlapMatrix(annoM$utr5,peak)

  utr5Mat <- Reduce("+",lapply(dataUTR5,calCovWin,annoS=annoM$utr5S,anno=annoM$utr5,mQuery=annoM$utr5M,fra=frac[1]))
  cdsMat <- Reduce("+",lapply(dataCDS,calCovWin,annoS=annoM$cdsS,anno=annoM$cds,mQuery=annoM$cdsM,fra=frac[2]))
  utr3Mat <- Reduce("+",lapply(dataUTR3,calCovWin,annoS=annoM$utr3S,anno=annoM$utr3,mQuery=annoM$utr3M,fra=frac[3]))

  mat <- as.data.frame(cbind(c(1:sum(frac)),c(utr5Mat,cdsMat,utr3Mat),c(rep("UTR5",frac[1]),rep("CDS",frac[2]),rep("UTR3",frac[3]))))

  colnames(mat) <- c("win","num","group")
  mat$win <- as.numeric(as.character(mat$win))
  mat$num <- as.numeric(as.character(mat$num))
  mat$num <- mat$num/sum(mat$num)
  return(mat)
}

##
calCovWin <- function(data,anno=cds,annoS=cdsS,mQuery=cdsM,fra=frac[2])
{
  s1 <- mQuery[[data[["subjectHits"]]]]
  wid <- max(max(data[["startP"]],data[["start"]])-data[["start"]],0)+1

  if( data[["queryHits"]] - s1 > 0 ) {
    tmp <- anno[subjectHits==data[["subjectHits"]]]
    left <- sum(tmp[tmp$queryHits<data[["queryHits"]],4]) + wid
  } else {
    left <- wid
  }
  right <- left+min(data[["end"]],data[["endP"]])-max(data[["start"]],data[["startP"]])
  all.wid <- annoS[[data[["subjectHits"]]]]

  left.bin1 <- ceiling(fra*left/all.wid)
  right.bin1 <- ceiling(fra*right/all.wid)
  left.bin2 <- ifelse(data[["strand"]]=="+",left.bin1,(fra+1-right.bin1))
  right.bin2 <- ifelse(data[["strand"]]=="+",right.bin1,(fra+1-left.bin1))

  if(right.bin2==left.bin2) {
    winCounts = as.vector(rep(0,fra))
    weight <- (right-left)/(all.wid/fra)
    winCounts[right.bin2]=winCounts[right.bin2]+ weight
  } else {
    winCounts = as.vector(rep(0,fra))
    weight.left <- left.bin1-left/(all.wid/fra)
    weight.right <- right/(all.wid/fra)-(right.bin1-1)
    winCounts[left.bin2]=winCounts[left.bin2]+weight.left
    winCounts[right.bin2]=winCounts[right.bin2]+weight.right
  }

  if(right.bin2-left.bin2 > 1)  {
    winCounts[(left.bin2+1):(right.bin2-1)]=winCounts[(left.bin2+1):(right.bin2-1)]+1
  }
  return(winCounts)
}

##

orderAndCom <- function(anno,longestIsoform) {
  anno <- as.data.frame(reduce(anno))
  anno <-anno[order(anno$seqnames,anno$start),]
  ans <- findOverlaps(GRanges(anno),GRanges(longestIsoform))
  loci <- subsetByOverlaps(GRanges(anno),GRanges(longestIsoform))

  res <- cbind(as.data.frame(loci),as.data.frame(ans))
  res$subjectHits<-factor(res$subjectHits)
  levels(res$subjectHits) <- 1:length(levels(res$subjectHits))
  res <- as.data.table(res)
  res <- res[order(res$subjectHits,res$queryHits),]
  return(res)
}

##

overlapMatrix <- function(anno.new,peak) {
  ov <- suppressWarnings(findOverlaps(GRanges(anno.new),GRanges(peak)))
  p <- Pairs(GRanges(anno.new), GRanges(peak), hits=ov)
  p <- as.data.frame(p)[,c(1:7,11:13)]
  colnames(p) <- c("seqnames","start","end","width","strand","queryHits","subjectHits","startP","endP","widthP")
  ll <- vector("list",nrow(p))
  for (i in 1:nrow(p)) { ll[[i]] <- p[i,] }
  return(ll)
}

