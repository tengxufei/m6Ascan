#' @title glmTest
#' @description significant test
#' @param counts counts file generated from previous step
#' @param condition The input bam file.
#' @param cutoff The gtf format gene annotation file
#' @param rep This option is the same in featureCounts, if your data is strandSpecfic, this value should be 1(stranded) or 2(reversely stranded), if not, it should be 0.
#' @import data.table
#' @export


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



