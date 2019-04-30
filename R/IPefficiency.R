# Most scripts here is copied from Trumpet:https://github.com/skyhorsetomoon/Trumpet
#' @title IPefficiency
#' @param IP.bam The ip bam file.
#' @param Input.bam The input bam file.
#' @param gtf The gtf format gene annotation file
#' @param sampleName your sample name, only used for names of print plot and table
#' @import GenomicFeatures
#' @export

effiTest <- function(IP.bam=c(f2,f4,f6),
                     Input.bam=c(f1,f3,f5),
                     gtf=NULL,
                     sampleName=NULL)
{

  result <- get_readscount(IP_BAM = IP.bam, Input_BAM = Input.bam,GENE_ANNO_GTF = gtf)

  s <- result[[1]]; ind <- unique(s$pos); len <- length(ind); n <- nrow(s)
  se <- seq(1, n, len); sa <- s[, -(1:2)]
  sample_name <- get.sampleid(IP_BAM = ip_bam, Input_BAM= input_bam)
  IP_groupname <- sample_name[[1]]
  Input_groupname <- sample_name[[2]]

  group_IP <- as.matrix(sa[, (seq_len(length(IP_groupname)))])
  Group_Input <- sa[, -(seq_len(length(IP_groupname)))]
  if(is.null(ncol(Group_Input))) {
    group_Input <- as.matrix(Group_Input)
  } else {
    group_Input <- Group_Input[, -((length(Input_groupname) + 1):ncol(Group_Input))]
    group_Input <- as.matrix(group_Input)
  }
  Unit_Input <- Unit_bam(group_Input,se,ind,len,n)
  Unit_IP <- Unit_bam(group_IP,se,ind,len,n)

  out <- .SES_IP(group_IP, group_Input, IP_groupname,se,len,n)
  newa <- out[[1]]; a <- out[[2]]; Enrich_table <- out[[3]]
  vline1 <- data.frame(ID = IP_groupname, pos = a)

  colnames(newa)<-c("pos","Sample","pro","ID")
  pos <- newa$pos; pro <- newa$pro; Sample <- newa$Sample
  pdf(paste(sampleName,"-efficiency-evaluation.pdf"),width=7,height=4)
  ggplot(data = newa, aes(x = pos, y = pro, colour = Sample )) +
    geom_line() + facet_grid(ID~.) +
    geom_vline(aes(xintercept = pos), vline1) + theme_bw()+
    theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
          panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"),axis.text = element_text(size = 12,color = "black"),
          title = element_text(size = 9),
          plot.title = element_text(hjust = 0.5),
          legend.key.height=unit(0.5,'cm'),
          legend.key.width=unit(0.25,'cm'),
          legend.text=element_text(size=9),
          legend.title=element_text(size=9))+
    labs(x = "Percentage of Bins", y = "Percentage of Reads")
  dev.off()

  colnames(Enrich_table) <- c("Sample ID", "Percent of Region Enriched with Signal", "Scale Factor")
  write.table(Enrich_table,paste(sampleName,"-efficiency-evaluation.xls"),sep = "\t",row.names = F)

}

.normalize_sample<-function(s1){
  row.sum<-rowSums(s1)
  z<-which(row.sum>10)
  s2<-s1[z,]
  row.mean<-rowMeans(s2)
  for(i in seq_len(length(row.mean))){
    if(row.mean[i]<2)
      row.mean[i]<-2
  }
  s3<-apply(s2,2,function(x,a)x/a, a=row.mean)
  s3<-as.matrix(s3)
  return(s3)
}

.singleBAMreads <- function(bam, se, len, n) {
  w <- vector(mode = "numeric", length = 0)
  for (i in seq_len(len)) {
    se <- seq(i, n, len)
    w <- cbind(w, bam[se])
  }
  return(w)
}

.unified_sample<-function(group,se,ind,len,n)
{
  m<-matrix(nrow = length(se),ncol = length(ind))
  v<-matrix(data=0,nrow=length(se),ncol = length(ind))
  group<-as.matrix(group)
  for(i in seq_len(ncol(group))){
    m<-.singleBAMreads(group[,i],se,len,n)
    v<-m+v
  }
  v<-v/length(ncol(group))
  return(v)
}

.Norm_bam <- function(bam, se, len, n) {
  s1 <- .singleBAMreads(bam, se, len, n)
  s3 <- .normalize_sample(s1)
  bam <- sapply(t(s3), unlist)
  return(bam)
}  ##get IP


Unit_bam <- function(group,se,ind,len,n) {
  v <- .unified_sample(group,se,ind,len,n)
  Input1 <- .normalize_sample(v)
  Input <- sapply(t(Input1), unlist)
  return(Input)
}

.SES_IP <- function(group_bam, unit_bam, IP_group_name, se, len, n) {
  new <- data.frame()
  a <- b <- z <- vector(mode = "numeric", length = 0)
  ID <- vector()
  for (i in seq_len(length(IP_group_name))) {

    bam <- vector(mode = "numeric", length = 0); bam <- .Norm_bam(group_bam[, i], se, len, n)
    Input1 <- unit_bam
    if(length(bam)==length(Input1)){
      ip1 <- bam; InPut1 <- Input1
    }
    M <- max(length(bam), length(Input1))
    if(length(bam) < M && length(Input1)==M){
      bam_sample <- sample(bam, size = (M-length(bam)),replace = TRUE)
      ip1 <- c(bam_sample, bam)
      InPut1 <-  Input1
    }
    if(length(Input1) < M && length(bam)==M){
      Input1_sample <- sample(Input1, size = (M-length(Input1)),replace = TRUE)
      ip1 <- bam
      InPut1 <- c(Input1_sample, Input1)
    }
    v1 <- sort(ip1); v2 <- sort(InPut1)
    x <- v1 - min(v1); x1 <- v2 - min(v2)
    ip <- x/sum(x); Input <- x1/sum(x1)
    cum_bam <- vector(mode = "numeric", length = 0); cum_bam <- cumsum(ip)
    newpos <- 1:length(ip); pos <- newpos/length(ip)
    unified_Input <- vector(mode = "numeric", length = 0); unified_Input <- cumsum(Input)
    c <- (unified_Input - cum_bam)
    a[i] <- pos[which.max(c)]; z[i] <- max(c)
    com <- as.data.frame(cbind(pos, cum_bam, unified_Input))
    new1 <- as.data.frame(melt(data = com, id = "pos", value.name = "pro"))
    var <- rep(c(IP_group_name[i], "unified_Input"), c(nrow(new1[(new1$variable) ==
                                                                   "cum_bam", ]), nrow(new1[(new1$variable) == "unified_Input", ])))
    new1$variable <- var
    b[i] <- length(new1$pos)
    new <- rbind(new, new1, data.frame())
    ID1 <- rep(IP_group_name[i], b[i])
    ID <- c(ID, ID2, vector())
  }
  new <- as.data.frame(cbind(new, ID))
  Scale_factor <- round(z, 2)
  p <- vector()
  for (i in seq_len(length(IP_group_name))) {
    p[i] <- paste(round((1 - a[i]) * 100, 2), "%")
  }
  Sample <- IP_group_name
  Enrichment_region <- p
  Enrich_table <- as.data.frame(cbind(Sample, Enrichment_region, Scale_factor))
  return(list(new, a, Enrich_table))
}

get_readscount <- function(IP_BAM, Input_BAM,
                           GENE_ANNO_GTF = NA, GENOME = NA, UCSC_TABLE_NAME = "knownGene", TXDB = NA,
                           sample_size = NA) {

  # download the annotation
  if (suppressWarnings((!is.na(GENOME)) & (!is.na(UCSC_TABLE_NAME)) &
                       is.na(TXDB) & is.na(GENE_ANNO_GTF))) {
    op <- options(warn = (-1))
    txdb = makeTxDbFromUCSC(genome = GENOME, tablename = UCSC_TABLE_NAME)
    options(op)
  }
  if (suppressWarnings(!is.na(GENE_ANNO_GTF) & is.na(TXDB))) {
    op <- options(warn = (-1))
    txdb <- makeTxDbFromGFF(GENE_ANNO_GTF, format = "gtf")
    options(op)
  }

  # use provided annotation data file
  if (suppressWarnings(!is.na(TXDB))) {
    txdb <- loadDb(TXDB)
  }
  gc <- makeGuitarCoordsFromTxDb(txdb, noBins = 20)
  gc_info <- mcols(gc)
  file <-  index <- as.character(c(IP_BAM, Input_BAM))
  sample_name <- as.vector(unlist(get.sampleid(IP_BAM, Input_BAM)))
  transform_table <- as.data.frame(cbind(file, sample_name))
  colnames(transform_table) <- c("files", "sample ID")
  # get reads count
  result2 <- gc_info
  noFiles <- length(file)
  total_reads <- vector(length = noFiles)
  exon_reads <- vector(length = noFiles)
  intron_reads <- vector(length = noFiles)
  no_genic <- vector(length = noFiles)
  percent_exon <- vector(length = noFiles)
  percent_intron <- vector(length = noFiles)
  percent_nogenic <- vector(length = noFiles)
  UTR5_reads <- vector(length = noFiles)
  CDS_reads <- vector(length = noFiles)
  UTR3_reads <- vector(length = noFiles)
  percent_UTR5 <- vector(length = noFiles)
  percent_CDS <- vector(length = noFiles)
  percent_UTR3 <- vector(length = noFiles)
  if (is.na(sample_size)) {
    for (i in seq_len(noFiles)) {
      print(paste("working on the ", i, "-th bam file ...", sep = ""))
      bam <- readGAlignments(file[i])
      total_reads[i] <- paste0(round(length(bam)/10^6, 2), "M")
      bin_count <- data.frame(countOverlaps(gc, bam))
      names(bin_count) <- sample_name[i]
      result2 <- data.frame(result2, bin_count)
      exon <- exonsBy(txdb, by = "tx")
      exon_count <- countOverlaps(bam, exon)
      exon_reads[i] <- paste0(round(sum(exon_count>0)/10^6, 2),
                              "M")
      percent_exon[i] <- paste0(round((sum(exon_count>0)/(length(bam)))*100,2), "%")
      percent_exon[i] <- paste0("(", percent_exon[i], ")")

      intron <- intronsByTranscript(txdb)
      intron_count <- countOverlaps(bam, intron)
      intron_reads[i] <- paste0(round(sum(intron_count > 0)/10^6,
                                      2), "M")
      percent_intron[i] <- paste0(round((sum(intron_count > 0)/(length(bam)))*100,2), "%")
      percent_intron[i] <- paste0("(", percent_intron[i], ")")

      no_genic[i] <- paste0(round((length(bam)-sum(exon_count>0)-sum(intron_count>0))/10^6,2), "M")
      percent_nogenic[i] <- paste0(round(((length(bam)-sum(exon_count>0)-sum(intron_count>0))/length(bam))*100,2), "%")
      percent_nogenic[i] <- paste0("(", percent_nogenic[i], ")")

      utr5 <- fiveUTRsByTranscript(txdb)
      utr5_count <- countOverlaps(bam, utr5)
      UTR5_reads[i] <- paste0(round(sum(utr5_count > 0)/10^6, 2),
                              "M")
      cds <- cdsBy(txdb, by = "tx")
      cds_count <- countOverlaps(bam, cds)
      CDS_reads[i] <- paste0(round(sum(cds_count > 0)/10^6, 2), "M")
      utr3 <- threeUTRsByTranscript(txdb)
      utr3_count <- countOverlaps(bam, utr3)
      UTR3_reads[i] <- paste0(round(sum(utr3_count > 0)/10^6, 2),
                              "M")
      sum_component <- sum(utr5_count > 0) + sum(cds_count > 0) +
        sum(utr3_count > 0)
      percent_UTR5[i] <- paste0(round((sum(utr5_count > 0)/(sum_component)) *
                                        100, 2), "%")
      percent_UTR5[i] <- paste0("(", percent_UTR5[i], ")")
      percent_CDS[i] <- paste0(round((sum(cds_count > 0)/(sum_component)) *
                                       100, 2), "%")
      percent_CDS[i] <- paste0("(", percent_CDS[i], ")")
      percent_UTR3[i] <- paste0(round((100 - round((sum(utr5_count > 0)/(sum_component)) *
                                                     100, 2) - round((sum(cds_count > 0)/(sum_component)) *
                                                                       100, 2)),2), "%")
      percent_UTR3[i] <- paste0("(", percent_UTR3[i], ")")
    }
  }
  if (!is.na(sample_size)) {
    sample_size <- as.numeric(sample_size)
    for (i in seq_len(noFiles)) {
      print(paste("working on the ", i, "-th bam file ...", sep = ""))
      bam <- readGAlignments(file[i])
      noR <- length(bam)
      id <- sample.int(noR, size = sample_size, replace = TRUE)
      bam <- bam[id]
      total_reads[i] <- paste0(round(length(bam)/10^6, 2), "M")
      bin_count <- countOverlaps(gc, bam)
      bin_count <- data.frame(bin_count)
      names(bin_count) <- sample_name[i]
      result2 <- data.frame(result2, bin_count)
      exon <- exonsBy(txdb, by = "tx")
      exon_count <- countOverlaps(bam, exon)
      exon_reads[i] <- paste0(round(sum(exon_count > 0)/10^6, 2),
                              "M")

      percent_exon[i] <- paste0(round((sum(exon_count>0)/(length(bam)))*100,2), "%")
      percent_exon[i] <- paste0("(", percent_exon[i], ")")

      intron <- intronsByTranscript(txdb)
      intron_count <- countOverlaps(bam, intron)
      intron_reads[i] <- paste0(round(sum(intron_count > 0)/10^6,
                                      2), "M")
      percent_intron[i] <- paste0(round((sum(intron_count > 0)/(length(bam)))*100,2), "%")
      percent_intron[i] <- paste0("(", percent_intron[i], ")")

      no_genic[i] <- paste0(round((length(bam)-sum(exon_count>0)-sum(intron_count>0))/10^7,2), "M")
      percent_nogenic[i] <- paste0(round(((length(bam)-sum(exon_count>0)-sum(intron_count>0))/length(bam))*100,2), "%")
      percent_nogenic[i] <- paste0("(", percent_nogenic[i], ")")

      utr5 <- fiveUTRsByTranscript(txdb)
      utr5_count <- countOverlaps(bam, utr5)
      UTR5_reads[i] <- paste0(round(sum(utr5_count > 0)/10^6, 2),
                              "M")
      cds <- cdsBy(txdb, by = "tx")
      cds_count <- countOverlaps(bam, cds)
      CDS_reads[i] <- paste0(round(sum(cds_count > 0)/10^6, 2), "M")
      utr3 <- threeUTRsByTranscript(txdb)
      utr3_count <- countOverlaps(bam, utr3)
      UTR3_reads[i] <- paste0(round(sum(utr3_count > 0)/10^6, 2),
                              "M")
      sum_component <- sum(utr5_count > 0) + sum(cds_count > 0) +
        sum(utr3_count > 0)
      percent_UTR5[i] <- paste0(round((sum(utr5_count > 0)/(sum_component)) *
                                        100, 2), "%")
      percent_UTR5[i] <- paste0("(", percent_UTR5[i], ")")
      percent_CDS[i] <- paste0(round((sum(cds_count > 0)/(sum_component)) *
                                       100, 2), "%")
      percent_CDS[i] <- paste0("(", percent_CDS[i], ")")
      percent_UTR3[i] <- paste0(round((100 - round((sum(utr5_count > 0)/(sum_component)) *
                                                     100, 2) - round((sum(cds_count > 0)/(sum_component)) *
                                                                       100, 2)),2), "%")
      percent_UTR3[i] <- paste0("(", percent_UTR3[i], ")")
    }

  }
  reads_exon <- paste(exon_reads, percent_exon)
  reads_intron <- paste(intron_reads, percent_intron)
  reads_nogenic <- paste(no_genic, percent_nogenic)
  reads_UTR5 <- paste(UTR5_reads, percent_UTR5)
  reads_CDS <- paste(CDS_reads, percent_CDS)
  reads_UTR3 <- paste(UTR3_reads, percent_UTR3)
  read_alignment_summary <- cbind(sample_name, total_reads,reads_exon, reads_intron, reads_nogenic, reads_UTR5, reads_CDS, reads_UTR3)
  read_alignment_summary <- as.data.frame(read_alignment_summary)
  colnames(read_alignment_summary) <- c("sample ID","total reads#", "exon reads", "intron reads", "no genic reads", "5'UTR reads", "CDS reads", "3'UTR reads")
  # remove DNA regions
  i <- which(result2$comp == "Front")  # remove front DNA
  result2 <- result2[-i, ]
  i <- which(result2$comp == "Back")  # remove tail DNA
  result <- result2[-i, ]
  t <- data.frame(result)
  t1 <- t[(t$category) == "mRNA", ]
  t2 <- t1[(t1$comp == "UTR5"), ]
  t3 <- t1[(t1$comp == "CDS"), ]
  t3$pos <- t3$pos + 1
  t4 <- t1[(t1$comp == "UTR3"), ]
  t4$pos <- t4$pos + 2
  t0 <- rbind(t2, t3, t4)
  s <- data.frame()
  s <- aggregate(cbind(t0[, 6]) ~ pos + txid, t0, mean)
  for (i in (length(t0) - noFiles + 2):length(t0)) {

    w <- aggregate(cbind(t0[, i]) ~ pos + txid, t0, mean)
    w <- w[, -(1:2)]
    s <- cbind(s, w)
  }
  s <- as.data.frame(s)
  colnames(s) <- c("pos", "txid", sample_name)
  read_result <- list(s, read_alignment_summary, transform_table)
  return(read_result)
}

get.sampleid <- function(IP_BAM, Input_BAM) {

  IP_name <- paste0("IP", seq_along(IP_BAM))
  Input_name <- paste0("Input", seq_along(Input_BAM))

  sample_name <- list()
  sample_name$IP_name <- IP_name
  sample_name$Input_name <- Input_name
  return(sample_name)
}

