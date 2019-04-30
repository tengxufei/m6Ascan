#' @title normalizeLibrary
#' @param counts counts file generated from previous step
#' @param input.site The input sites in counts dataframe.
#' @param IP.site The IP sites in counts dataframe.
#' @import DESeq2
#' @export


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



