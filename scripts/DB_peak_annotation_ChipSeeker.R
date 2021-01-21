library(ChIPseeker)
library(ReactomePA)
library(clusterProfiler)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# annotation db - use 5kb window
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGenes
promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000)

# IDR merging of peaks called with MACS2 p < 0.05
# split human peaks from merged MACS2 narrow.peaks output before IDR call

# read MYCN (m) and MYCN + DOX (d) peaks
dpeaks <- readPeakFile('~/scratch/chipseq_wf/all_bam_p05/human/idr_MYCN_DOX.peaks')
mpeaks <- readPeakFile('~/scratch/chipseq_wf/all_bam_p05/human/idr_MYCN.peaks')


# blacklisted files
dpeaks <- readPeakFile('~/scratch/chipseq_wf/peaks/human/idr_MYCN_DOX.peaks')
#dpeaks <- dpeaks[dpeaks$X1000 >= 540, ]
mpeaks <- readPeakFile('~/scratch/chipseq_wf/peaks/human/idr_MYCN.peaks')
#mpeaks <- mpeaks[mpeaks$X1000 >= 540, ]

# get matrix / counts for plot
dtagMatrix <- getTagMatrix(dpeaks, windows=promoter)
mtagMatrix <- getTagMatrix(mpeaks, windows=promoter)

dtagcount <- ChIPseeker:::getTagCount(dtagMatrix, xlim=c(-5000, 5000))
mtagcount <- ChIPseeker:::getTagCount(mtagMatrix, xlim=c(-5000, 5000))



#------- plot
par(mar=c(6,8,2,2))
plot(dtagcount, type='l', 
     col=rgb(0, 0, 1, 0.5), 
     lwd=2, bty='n', las=1,
     xlab="Genomic Region (5'->3')", 
     ylab = "Read Count Frequency", 
     mgp=c(4.5,1,0))
lines(mtagcount, col=rgb(1, 0, 0, 0.3), lwd=1.5)
abline(v=0, col=rgb(0.3, 0.3, 0.3, 0.3))
legend(2000, 0.0002, pch=15, c('MYCN', 'MYCN + DOX'), 
       cex=0.8, bty='n',
       col=c(rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.5)))

dev.copy2pdf(file='~/scratch/chipseq_wf/Chipseq_mycn/figures/binding_relative_to_TSS_5KB.pdf', height=3, width=6)

#tagHeatmap(dtagMatrix)

# annotations
dpeakAnno <- annotatePeak(dpeaks, tssRegion=c(-5000, 5000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")
mpeakAnno <- annotatePeak(mpeaks, tssRegion=c(-5000, 5000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")
# get 
dgids <- unique(as.data.frame(dpeakAnno)$SYMBOL)
dim(data.frame(dgids))

old_gids <- read.csv('~/scratch/chipseq_wf/Chipseq_mycn/data/IDR_DOX_peak_set_annotations.csv')

sum(dgids %in% old_gids$x)

write.csv(dgids, row.names = FALSE, file = '~/scratch/chipseq_wf/Chipseq_mycn/data/IDR_DOX_peak_set_annotations_BL.csv', quote=FALSE)
mmgids <- unique(as.data.frame(mpeakAnno)$SYMBOL)



length(mgids)
length(dgids)

diffentrez <- setdiff(as.data.frame(dpeakAnno)$geneId, as.data.frame(mpeakAnno)$geneId)
diffgids <- setdiff(as.data.frame(dpeakAnno)$SYMBOL, as.data.frame(mpeakAnno)$SYMBOL)
grep('STAT', dgids)




sum(dgids %in% mgids)
sum(mgids %in% dgids)

dpathway <- enrichPathway(as.data.frame(dpeakAnno)$geneId, readable = TRUE)
mpathway <- enrichPathway(as.data.frame(mpeakAnno)$geneId, readable = TRUE)
diffpathway <- enrichPathway(diffentrez, readable = TRUE)

dres <- dpathway@result[dpathway@result$p.adj < 0.05, ]
mres <- mpathway@result[mpathway@result$p.adj < 0.05, ]

danno <- as.data.frame(dpeakAnno)$geneId

ego <- enrichGO(gene = danno, 
                keyType = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

ego <- enricher(danno)

res_ego <- ego@result[ego@result$p.adjust < 0.05, ]

res_ego$geneID[res_ego$Description=='regulation of toll-like receptor signaling pathway']

#---------

degs <- read_excel('~/MYCN DEG Limma analysis.xlsx')
head(degs)

sum(degs$`MYCN downregulated genes` %in% dgids)
sum(degs$`MYCN upregulated genes` %in% dgids)
sum(degs$`MYCN downregulated genes` %in% mgids)
sum(degs$`MYCN upregulated genes` %in% mgids)

#------------ motifs

library(motifmatchr)
library(GenomicRanges)
library(SummarizedExperiment)

# load some example motifs
data(example_motifs, package = "motifmatchr") 
library(BCRANK)
library(Biostrings)
library(seqLogo)

fpeaks <- fread('~/scratch/chipseq_wf/all_bam_p05/human/idr_MYCN_DOX_all.peaks')
fpeaks <- fpeaks[rev(order(V12)), ]

pk_ord <- paste0(fpeaks$V1, ':', fpeaks$V2, '-', fpeaks$V3)

head(fpeaks)

pk_fa <- readDNAStringSet('~/scratch/chipseq_wf/motifs/MYN_DOX_peak_seqs.fa')
pk_fa_ord <- pk_fa[pk_ord]
bcrank_out <- bcrank('~/scratch/chipseq_wf/motifs/MYN_DOX_peak_seqs.fa', restarts=25, use.P1 = TRUE, use.P2 = TRUE)

identical(names(pk_fa_ord), pk_ord)

writeXStringSet(pk_fa_ord, '~/scratch/chipseq_wf/motifs/MYCN_DOX_sorted_peak_seqs.fa')
bcrank_out <- bcrank('~/scratch/chipseq_wf/motifs/MYCN_DOX_sorted_peak_seqs.fa', restarts=20, use.P1 = TRUE, use.P2 = TRUE, plot.progress = TRUE)

bcrank_out <- bcrank('~/scratch/chipseq_wf/motifs/MYCN_DOX_peak_seqs_all.fa', restarts=25, use.P1 = TRUE, use.P2 = TRUE, plot.progress = TRUE)


topMotif <- toptable(bcrank_out, 2)
weightMatrix <- pwm(topMotif, normalize = FALSE)
weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
seqLogo(weightMatrixNormalized)

hist(fpeaks$V5)


head(names(pk_fa))


sum(!is.na(degs$`MYCN upregulated genes`))

455/18000

47/455
1576/18000


length(degs$`MYCN downregulated genes`)
