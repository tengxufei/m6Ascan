#' @title segment
#' @description segmentation algorithm and then map transcript loci from genome loci
#' @param gtf The splited gtf file
#' @param catMethyl The output file of sig.test()
#' @param threads Threads used for parallel running
#' @import data.table
#' @export


reportResults <- function(gtf=win.gtf, catMethyl=NULL, threads=4)
{
  gtf <- fread(gtf)[,c(1,4,5,7,9)]
  colnames(gtf) <- c("chr","start","end","strand","geneBin")
  x <- data.table(gtf[,-5],"ID"=apply(gtf[,5],1,function(x){gsub("\";","",gsub("bin_id \"","",x))}))
  new.gtf <- separate(data = x, col = ID, into = c("gene", "bin"), sep = "-")

  catMethyl <- fread(catMethyl)

  geneName <- levels(factor(catMethyl$gene)); no.genes=length(geneName)

  registerDoParallel( cores = threads)

  ## narrow peaks and map to transcript loci
  narrowAndMap <- foreach(i = 1:no.genes) %dopar% {

    df <- catMethyl[catMethyl$gene==geneName[i],]
    colN <- (ncol(df)-4)/3
    input.site <- c(3:(3+2*colN-1))[c(3:(3+2*colN-1))%% 2 != 0]

    interval <- mean(df$bin[2:length(df$bin)]- df$bin[1:(length(df$bin)-1)])
    binNum <- ifelse(interval < 2, 10, ifelse(interval >= 5, 3, 5))
    # df$IP1 or mean(df$IP1,df$IP2)
    big <- findMaxMean(df$IP1,binNum,nrow(df))
    df <- df[big:(big+binNum-1),]

    tmp.gtf <- new.gtf[(new.gtf$gene==geneName[i]) && (new.gtf$bin %in% df$bin),]
    ## map transcript loci
    block <- mapToTrans(tmp.gtf,step=step)
    maxRow <- which.max(c(rowMeans(df[,(input.site+1)])))

    data.table("chrom" = tmp.gtf[1,1],
               "chromStart" = min(tmp.gtf[tmp.gtf$bin==df[1,2],"start"]),
               "chromEnd" = max(tmp.gtf[tmp.gtf$bin==df[nrow(df),2],"end"]),
               "name" = geneName[i],
               "score" = ".",
               "strand" = tmp.gtf[1,4],
               "thickStart" = min(tmp.gtf[tmp.gtf$bin==df[1,2],"start"]),
               "thickend" = max(tmp.gtf[tmp.gtf$bin==df[nrow(df),2],"end"]),
               "itemRgb" = 0,
               "blockCount" = block[[1]],
               "blockSizes" = block[[2]],
               "blockStarts" = block[[3]],
               "mean.Input.value" = mean(c(df[,input.site])),
               "mean.IP.value" = mean(c(df[,(input.site+1)])),
               "mean.max.Input.value" = mean(df[maxRow,input.site]),
               "mean.max.IP.value" = mean(df[maxRow,input.site]),
               "mean.max.ES" = mean(df[maxRow,c(-1:-colN)]),
               "max.pVals" = df[maxRow,"pVals"],
               "max.fdr" = df[maxRow,"edgeR_fdrs"]
    )
  }

  return(rbindlist(narrowAndMap))

}


mapToTrans <- function(tmp.gtf, step=10)
{
  y <- tmp.gtf[1,1]; row <- 1
  for (i in 1:(nrow(tmp.gtf)-1)) {
    if((tmp.gtf[i+1,2]-tmp.gtf[i,2]) > step){y <- c(y,tmp.gtf[i+1,2]);row=c(row,i+1)}
  }

  if(length(y)==1) {size=sum(tmp.gtf[,2]-tmp.gtf[,1])}
  else {
    row2 <- c((row[-1]-1),nrow(tmp.gtf))
    a2 <- as.data.frame(cbind(row,row2)) %>% rowwise() %>% mutate(size = sum(tmp.gtf[row:row2,2]-tmp.gtf[row:row2,1]))[,"size"]
    size <- a2$size
  }

  return(list(length(y),size,y))
}



