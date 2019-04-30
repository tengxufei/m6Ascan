#' @title PeakCalling
#' @param IP.bam The ip bam file.
#' @param Input.bam The input bam file.
#' @param gtf The gtf format gene annotation file
#' @param strandSpecific This option is the same in featureCounts, if your data is strandSpecfic, this value should be 1(stranded) or 2(reversely stranded), if not, it should be 0.
#' @param isPairedEnd logical. Depends on your bam file is pair end or single end.
#' @param threads Threads used for parallel running
#' @import Rsubread
#' @import edgeR
#' @import tidyr
#' @import dplyr
#' @import statmod
#' @import magrittr
#' @export

get.win.counts <-function(IP.bam=c(f2,f4,f6,f8),
                          Input.bam=c(f1,f3,f5,f7),
                          gtf=win.gtf,
                          strandSpecific=0,
                          isPairedEnd=F,
                          threads=4
)
{

  winCounts <- featureCounts(files=cat.sample(Input.bam,IP.bam),annot.ext=gtf,
                             isGTFAnnotationFile=TRUE, GTF.featureType="bin", GTF.attrType="bin_id",
                             isPairedEnd=isPairedEnd, strandSpecific=strandSpecific, allowMultiOverlap=T,
                             countMultiMappingReads=F,nthreads=threads)

  x = data.frame(winCounts$annotation[,c("GeneID")],winCounts$counts,stringsAsFactors=FALSE)
  colnames(x) <- c("ID",paste(rep(c("Input","IP"),(ncol(x)-1)/2),rep(c(1:((ncol(x)-1)/2)),each=2),sep=""))
  counts <- separate(data = x, col = ID, into = c("gene", "bin"), sep = "-")

  return(counts)
}

calc.NormFactors <- function(counts=win.counts, input.site=input.site, IP.site=IP.site) {

  input.median.mat <- sapply(counts[,input.site],function(x){tapply(x,counts$gene,median)})
  input.normF <- estimateSizeFactorsForMatrix(input.median.mat)

  enrich.median.mat <- sapply(counts[,IP.site],function(x){tapply(x,counts$gene,median)})/input.median.mat
  IP.normF <- estimateSizeFactorsForMatrix(enrich.median.mat)

  norm.factors <- vector()
  for (i in 1:length(IP.normF)) {
    norm.factors <- c(norm.factors,input.normF[i],IP.normF[i])
  }

  return(norm.factors)
}

cat.sample <- function(x1,x2) {
  b <- c()
  for (i in 1:length(x1)) {
    b <- c(b,x1[i],x2[i])
  }
  return(b)
}

calcES <- function(input.site, counts)
{
  medianRatio <- vector("list",length(input.site))
  for(i in input.site) {
    medianRatio[[(i-1)/2]] <- tapply(counts[,i],counts$gene,median)/tapply(counts[,i+1],counts$gene,median)
  }

  for(i in input.site) {
    counts[,paste("ES",(i-1)/2,sep="")] <- medianRatio[[(i-1)/2]][counts$gene]*counts[,i+1]/counts[,i]
  }
  return(counts)
}

# glm.nb test
glmTest <- function(counts) {
  group_list <- rep(c("Input","IP"),ncol(counts)/2-1)
  group_edgeR <- factor(group_list)
  # default design do not have other coefficient
  design <- model.matrix( ~ group_edgeR)
  # got the differential windows
  dge <- DGEList(counts=counts[,-1:-2], norm.factors =rep(1,ncol(counts)-2), group=group_list)
  # set the library size to 1 among all samples
  dge$samples$lib.size <- rep(1,ncol(counts)-2)
  # estimate dispersion based on default edgeR setting
  dge <- estimateDisp(dge, design = design, trend.method="none")
  # do the QLFit, if bad things happens, change to glmFit
  s <- try(fit <- glmQLFit(dge, design, robust=TRUE), silent=TRUE)
  if ('try-error' %in% class(s)) { fit <- glmFit(dge, design) }
  res <- glmLRT(fit)
  # get the pvalue and fdr
  pVals <- res$table[,4]
  names(pVals) <- rownames(res$table)
  edgeR_fdrs <- p.adjust(pVals, method = "fdr")
  # result
  glmGeneFit <- as.data.frame(data.table(counts,pVals,edgeR_fdrs))

  return(glmGeneFit)
}

# significant test
sig.test <- function(counts,
                     condition=c("KO","WT"),
                     cutoff=100,
                     rep=1)
{
  if(rep < 2) {
    cat("m6Ascan only used for data has biological replicates! Error is coming...\n")
  }
  n <- length(condition)

  for (i in 1:n) {
    tmp <- counts[,c(1:2,(3+(i-1)*rep*2):(3+i*rep*2-1),(3+n*rep*2+(i-1)*rep):(3+n*rep*2+i*rep-1))]
    tmp <- tmp[rowMeans(tmp[,c(4,6)]) > cutoff,]
    glmGeneFit <- data.table(glmTest(counts=tmp[,1:(rep+1)*2]),tmp[,((rep+1)*2+1):((rep+2)*2)])
    sigMethyl <- as.data.frame(glmGeneFit[glmGeneFit$edgeR_fdrs < 0.05 & glmGeneFit[,9] >2 & glmGeneFit[,10] >2, ])
    write.table(sigMethyl,file=paste(condition[i],"-catMethyl.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  }
}


