
library(edgeR)
library(EDASeq)
library(GenomicAlignments)
library(Hmisc)
library(dplyr)
library(ggplot2)
library(ChIPseeker)
library(txdbmaker)
library(org.Mm.eg.db)
library(ggupset)

## ---- paths  ----
counts_file <- "results_atacseq/bwa/merged_library/macs2/broad_peak/overlaps-1/counts_table/merged_peaks_macs2_overlap-1.counts.tsv"
peaks_bed   <- "results_atacseq/bwa/merged_library/macs2/broad_peak/overlaps-1/merged_overlap-1_peaksid.bed"
fasta       <- "references/Mus_musculus.GRCm38.dna_sm.toplevel.fa"
out_dir     <- "results/overlap-1/"
qc_dir      <- "results/overlap-1/QC/"
annot_base  <- "results/overlap-1/peak.annotation/"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)


######################################################################
##  Differential accessibility 
######################################################################

ff <- FaFile(fasta)

## load counts
cnt_table <- read.table(counts_file, sep = "\t", header = TRUE, blank.lines.skip = TRUE)
rownames(cnt_table) <- cnt_table$Geneid
colnames(cnt_table) <- gsub(".mLb.clN.sorted.bam", "", colnames(cnt_table))

# keep assembled chromosomes only
standard_chromosomes <- c(as.character(1:19), "X", "Y")
cnt_table <- cnt_table[cnt_table$Chr %in% standard_chromosomes, ]

# experimental groups (column order from featureCounts: HET then WT)
groups <- factor(c(rep("HET", 4), rep("WT", 3)))

# read counts to peaks (columns 7+)
reads.peak <- as.matrix(cnt_table[, 7:ncol(cnt_table)])

## GC content of peaks
gr <- GRanges(seqnames = cnt_table$Chr,
              ranges   = IRanges(cnt_table$Start, cnt_table$End),
              strand   = "*",
              mcols    = data.frame(peakID = cnt_table$Geneid))

# intersect with seqlevels(gr) so chromosomes absent from the data don't error
gr <- keepSeqlevels(gr, intersect(standard_chromosomes, seqlevels(gr)),
                    pruning.mode = "coarse")

peakSeqs       <- getSeq(ff, gr)
gcContentPeaks <- letterFrequency(peakSeqs, "GC", as.prob = TRUE)[, 1]
gcGroups       <- Hmisc::cut2(gcContentPeaks, g = 20)

# QC: count vs GC content per sample
lowList <- lapply(seq_len(ncol(reads.peak)), function(k) {
  set.seed(k); lowess(x = gcContentPeaks, y = log1p(reads.peak[, k]), f = 1/10)
})
dfAll <- do.call(rbind, lapply(seq_along(lowList), function(s) {
  o <- order(lowList[[s]]$x)
  data.frame(x = lowList[[s]]$x[o], y = lowList[[s]]$y[o], sample = colnames(reads.peak)[s])
}))
pdf(paste0(qc_dir, "GCcontent_peaks.pdf"))
print(ggplot(dfAll, aes(x, y, colour = sample)) + geom_line(linewidth = 1) +
        xlab("GC-content") + ylab("log(count + 1)") + theme_classic())
dev.off()

## GC-aware offsets
eda <- newSeqExpressionSet(
  counts      = reads.peak,
  featureData = AnnotatedDataFrame(
    data.frame(gc = gcContentPeaks, row.names = rownames(reads.peak))))
eda <- withinLaneNormalization(eda, "gc", num.bins = 20, which = "full", offset = TRUE)
eda <- betweenLaneNormalization(eda, which = "full", offset = TRUE)
dataOffset <- offst(eda)   

## edgeR with EDASeq offsets
groups <- relevel(groups, ref = "WT")
design <- model.matrix(~groups)

d    <- DGEList(counts = reads.peak, group = groups)
keep <- filterByExpr(d)
d    <- d[keep, , keep.lib.sizes = FALSE]

d$offset    <- -dataOffset[keep, ]
d.eda       <- estimateGLMCommonDisp(d, design = design)
fit         <- glmFit(d.eda, design = design)
lrt.EDASeq  <- glmLRT(fit, coef = 2)

DA_res <- as.data.frame(topTags(lrt.EDASeq, nrow(lrt.EDASeq$table)))
DA_res$Geneid <- rownames(DA_res)

## add coordinates and per-group median counts
count_df <- data.frame(
  Geneid            = rownames(d$counts),
  median_WT_counts  = apply(d$counts[, groups == "WT"], 1, median),
  median_HET_counts = apply(d$counts[, groups == "HET"], 1, median))

DA.res.coords <- DA_res %>%
  left_join(cnt_table[1:4], by = "Geneid") %>%
  left_join(count_df, by = "Geneid")

write.table(DA.res.coords, paste0(out_dir, "DA_HET_vs_WT.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

######################################################################
##  Peak annotation 
######################################################################

## Ensembl TxDb
txdb_ens <- makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL",
                                dataset = "mmusculus_gene_ensembl",
                                host    = "https://nov2020.archive.ensembl.org")


## helper: annotate a GRanges and write outputs
annotate_peaks <- function(gr, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ann <- annotatePeak(gr, tssRegion = c(-3000, 3000),
                      TxDb = txdb_ens, annoDb = "org.Mm.eg.db")

  df <- as.data.frame(ann)
  df$geneChr <- NULL
  colnames(df)[1:5] <- c("Chr", "peak_start", "peak_end", "peak_width", "peak_strand")
  write.table(df, file.path(out, "peaks_annotated.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(ann@annoStat, file.path(out, "annotation_summary.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  pdf(file.path(out, "AnnotVis.pdf")); print(upsetplot(ann, vennpie = TRUE)); dev.off()
  pdf(file.path(out, "TSSdist.pdf"), 6, 5)
  print(plotDistToTSS(ann, title = "ATAC-seq peaks relative to TSS")); dev.off()
  ann
}

## all consensus peaks
peaks.bed <- read.table(peaks_bed, sep = "\t", header = FALSE, blank.lines.skip = TRUE)
peaks.gr  <- GRanges(seqnames = peaks.bed[, 1],
                     ranges   = IRanges(peaks.bed[, 2], peaks.bed[, 3]),
                     strand   = "*",
                     mcols    = data.frame(peakID = peaks.bed[, 4]))

dir.create(annot_base, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(annot_base, "PeakCoverage.pdf"), 8, 20); covplot(peaks.gr); dev.off()
annotate_peaks(peaks.gr, annot_base)

## differentially accessible peaks (up / down)
to_gr <- function(x) GRanges(seqnames = x$Chr, ranges = IRanges(x$Start, x$End),
                             strand = "*", mcols = data.frame(peakID = x$Geneid))

annotate_peaks(to_gr(filter(DA.res.coords, FDR < 0.05, logFC > 0)),
               file.path(annot_base, "upregulated"))
annotate_peaks(to_gr(filter(DA.res.coords, FDR < 0.05, logFC < 0)),
               file.path(annot_base, "downregulated"))
