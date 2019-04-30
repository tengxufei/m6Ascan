#' @title plotSingle
#' @description visulization of MeRIP-seq data
#' @param site The bed file contains the region you want to plot.
#' @param gtf The gtf format gene annotation file
#' @param txdb Txdb file. Could be NULL when gtf is avaliable
#' @param annotation Annotation hub file
#' @param genome genome version, eg: "mm10", "hg38"
#' @param bwInfo File name of your *.bw file
#' @param bwDir The directory of your *.bw file
#' @import Gviz
#' @import RColorBrewer
#' @describeIn Attention: if you want to use txdb and annotation hub file, you need require the package first, eg: require(TxDb.Mmusculus.UCSC.mm10.knownGene); require(org.Mm.eg.db)
#' @export

plot.single.peak <- function(site,
                             gtf,
                             txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                             annotation=org.Mm.eg.db,
                             genome="mm10",
                             bwInfo,
                             bwDir="/your/work/directory"
)
{

  tryCatch(suppressPackageStartupMessages(library(Gviz)),
           error=(BiocManager::install("Gviz")))
  tryCatch(suppressPackageStartupMessages(library(RColorBrewer)),
           error=(install.packages("RColorBrewer")))

  if(is.character(gtf)) {
    cat("Making txdb from gtf file ...","\n")
    txdb=makeTxDbFromGFF(gtf,format="gtf")
  } else {
    cat("Getting txdb ...","\n")
    txdb = txdb
  }

grt <- GeneRegionTrack(txdb, genome=genome,showId=TRUE, geneSymbol=TRUE, name="UCSC")
z <- mapIds(annotation, gene(grt), "SYMBOL", "ENTREZID", multiVals = "first")

zz <- sapply(z, is.null); z[zz] <- gene(grt)[zz]
gr <- ranges(grt); mcols(gr)$symbol <- z
grt@range <- gr


# genome range
colnames(site)<-c("chr","start","end","strand")
chr<-site[rownames(site),]$chr
site$width<-with(site,end-start)
startpoint <- site[rownames(site),]$start
endpoint <- site[rownames(site),]$end


tracklist<-list()
# write chromosome to tracklist
itrack <- IdeogramTrack(genome = genome, chromosome = chr,outline=T)
tracklist[["itrack"]] <- itrack

# write scale bar to tracklist
scalebar <- GenomeAxisTrack(scale=0.25,col="black",fontcolor="black",name="Scale",labelPos="above",showTitle=TRUE)
tracklist[["scalebar"]] <- scalebar

#write genomic ranges to tracklist
axisTrack <- GenomeAxisTrack(labelPos="above",col="black",fontcolor="black",name=paste(chr,":",sep=""),exponent=0,showTitle=TRUE)
tracklist[["axisTrack"]] <- axisTrack

# write bigwig to tracklist
## color
colpal<-rep(brewer.pal(12,"Paired"),20)
coldf<-data.frame(col=colpal[1:nrow(bwInfo)],row.names = rownames(bwInfo),stringsAsFactors = F)

for(sample in rownames(bwInfo)){
  bgFile <- file.path(bwDir,paste(sample,".bw",sep=""))
  tracklist[[index]]<-DataTrack(range = bgFile,genome=genome,type="histogram", name=index,col.histogram=coldf[index,])
}

# write gene structure to tracklist
tracklist[["grt"]]<-grt


pdf("lociT.pdf",height=5,width=8)
plotTracks(tracklist, from = startpoint, to = endpoint,sizes = c(0.3,0.2,0.5,3/8,8/8,0.3),
           chromosome=chr,background.panel = "white", background.title = "white",
           col.title="black",col.axis="black",
           rot.title=0,cex.title=0.9,margin=38,title.width=1.5,
           collapseTranscripts = "longest")
dev.off()

}
