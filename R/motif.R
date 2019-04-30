#' @title segment
#' @description segmentation algorithm and then map transcript loci from genome loci
#' @describeIn  genome file should be required first, eg: require(BSgenome.Mmusculus.UCSC.mm10)
#' @param peakFile The m6A peak file
#' @param species The species, eg: "mm10", "hg38"
#' @param plotmethod The method to plot motif, "prob" or "bit"
#' @param sortCut Choose the most significant peaks to find motifs, this value means the amount of peaks you wanna choose
#' @import rGADEM
#' @import ggseqlogo
#' @export

motifFind <- function(peakFile=NULL,species=NULL,plotmethod="prob",sortCut=5000) {

  if(species=="Mmusculus") {
    require(BSgenome.Mmusculus.UCSC.mm10)
  }

  rgBED <- IRanges(start=peakFile[1:min(nrow(peakFile),sortCut),2],end=peakFile[1:min(nrow(peakFile),sortCut),3])
  Sequences<-RangedData(rgBED,space=peakFile[,1])
  gadem<-GADEM(Sequences,verbose=1,genome=species)
  motifList <- gadem@motifList
  resList <- list()
  for(i in 1:length(motifList)){
    resList[[i]] <- motifList[[i]]@pwm
  }

  ggseqlogo(resList[[1]], method=plotmethod)

}
