##Brundle workflow for sharing

library(Brundle)

##Initial conditions
jg.controlMinOverlap      <- 5
jg.treatedCondition       =  "DOX"
jg.untreatedCondition     =  "none"



###Read the files 

jg.ExperimentSampleSheet <- read.delim("~/path/to/github/Experiment_sample_sheet.csv", sep = ",",
                                       header = T)
jg.ControlSampleSheet <- read.delim("~/path/to/github/Control_sample_sheet.csv", sep = ",",
                                    header = T)

##Load the R objects generated in Brundle_step01.R


dbaExperiment <- readRDS("~/CHIP_seq/data_git/diffbind_exp")
dbaControl <- readRDS("~/CHIP_seq/data_git/diffbind_cont")

jg.experimentPeakset <- readRDS("~/CHIP_seq/data_git/control_counts_expt_object")
jg.controlPeakset <- readRDS("~/CHIP_seq/data_git/control_counts_control_object")

jg.controlCountsTreated <- readRDS("~/CHIP_seq/data_git/control_counts_treated_object")
jg.controlCountsUntreated <- readRDS("~/CHIP_seq/data_git/control_counts_untreated_object")

#Get the sample names for replicates that represent the two conditions.
jg.untreatedNames <- names(jg.controlCountsUntreated)
jg.treatedNames   <- names(jg.controlCountsTreated)

#Plot of the normalisation, for visualisation only, not necessary for analysis.

jg.plotNormalization(jg.controlCountsTreated,
                     jg.controlCountsUntreated)
## jg.plotNormalization finshes

#Calculate the noralisation coefficent
### Example of jg.getNormalizationCoefficient

jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                               jg.controlCountsUntreated)

#Calculates to correction factor for DiffBind
##Run by Finn
#jg.correctionFactor<-jg.getCorrectionFactor(jg.experimentSampleSheet,
#                                            jg.treatedNames,
#                                            jg.untreatedNames)
jg.correctionFactor <- readRDS("~/CHIP_seq/data_git/correctionFactor_object")

#Save data for examples in package
#save(jg.experimentPeakset, file="data/jg.experimentPeakset.rda")


#Apply the normalisation to the experimental peakset.

jg.experimentPeaksetNormalised<-jg.applyNormalisation(jg.experimentPeakset,
                                                      jg.coefficient,
                                                      jg.correctionFactor,
                                                      c("2a", "2b"))


#Return values to Diffbind and plot normalised result.
jg.dba <- DiffBind:::pv.resetCounts(dbaExperiment,
                                    jg.experimentPeaksetNormalised)

#make contrasts to compare DOX versus no DOX
jg.dba <- dba.contrast(jg.dba, jg.dba$masks$DOX, jg.dba$masks$none,
                       "DOX", "none")
jg.dba <-dba.analyze(jg.dba)
plot(jg.dba, contrast=1) ##sanity check for contrast and replicates


##Annotation of peak loci; This section can be tweaked further to alter the sensitivity of the analysis
library(ChIPpeakAnno)

# Loading TSS Annotation For Human Sapiens (GRCh38) Obtained From BiomaRt
data("TSS.human.GRCh38")

# Choosing the peaks for the interesting comparison, 
##create GRanges object from diffbind object
##Feel free to change fold to make the analysis more sensitive
jg.dba.gr <- dba.report(jg.dba, th=.05, bUsePval=TRUE, fold=0.5)
#data.peaks = dba.report(jg.dba, contrast=1)
#head(data.peaks)

# Annotate peaks with information on closest TSS using precompiled annotation data
data.peaksAnno=annotatePeakInBatch(jg.dba.gr, AnnotationData=TSS.human.GRCh38)
##use the following command to increase the span of annotation
#data.peaksAnno=annotatePeakInBatch(jg.dba.gr, AnnotationData=TSS.human.GRCh38, bindingRegion = c(-10000, 10000))

##convert ensembl to symbol and entrez (This conversion can be skipped if only Ensemble genes are to be used)
##see feature column of dataframe in line 106
library(biomaRt) 
ensembl_hs <- useMart(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      host = "www.ensembl.org")

# Add gene information
library(org.Hs.eg.db)
data.peaksAnno <- addGeneIDs(data.peaksAnno, org.Hs.eg.db, mart = ensembl_hs, feature_id_type = "ensembl_gene_id",
                             IDs2Add = c("symbol", "entrezgene"))
data.peaksAnno_df <- as.data.frame(data.peaksAnno)

mapping <- select(org.Hs.eg.db, 
                  keys = data.peaksAnno_df$symbol,
                  columns = c("ENTREZID", "SYMBOL"),
                  keytype = "SYMBOL")
data.peaksAnno_df$ENTREZID <- mapping[match(data.peaksAnno_df$symbol, mapping$SYMBOL),2]
write.table(data.peaksAnno_df, "~/CHIP_seq/data_git/DiffBind_Analysis_pval_th05_GRCh38.tsv", sep = "\t", quote = F, row.names = F)
##Annotation
#library(ChIPseeker)
#library(clusterProfiler)
library(ReactomePA)
pathways_peaks <- enrichPathway(data.peaksAnno_df$ENTREZID)
write.table(pathways_peaks, "~/CHIP_seq/data_git/Anno_DiffBind_Analysis_pval_th05_GRCh38.tsv", sep = "\t", quote = F, row.names = F)
dotplot(pathways_peaks)

#library(topGO)
library(enrichR)
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enricheddba <- enrichr(data.peaksAnno_df$symbol, dbs)
enricheddba_MF <- enricheddba[[1]]
enricheddba_CC <- enricheddba[[2]]
enricheddba_BP <- enricheddba[[3]]

enricheddba_BP_sig <- enricheddba_BP[enricheddba_BP$Adjusted.P.value < 0.05,]
write.table(enricheddba_BP_sig, "~/CHIP_seq/data_git/GO_BP_DiffBind_Analysis_adj_pval_th05_GRCh38.tsv", sep = "\t", quote = F, row.names = F)

##This concludes the DiffBind analysis

##################
##DEseq2 analysis and annotation
jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
jg.experimentPeaksetDeSeq<-jg.convertPeakset(jg.experimentPeakset)

jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)
jg.controlPeaksetDeSeq<-jg.convertPeakset(jg.controlPeakset)

### Establish size factors directly from Control data using DESeq2
jg.controlSizeFactors = estimateSizeFactorsForMatrix(jg.controlPeaksetDeSeq)

jg.conditions <- read.csv(file="~/CHIP_seq/data_git/Control_sample_sheet.csv", header=TRUE, sep=",")['Condition']
jg.controlDeSeq<-jg.runDeSeq(jg.controlPeaksetDeSeq, jg.conditions,jg.controlSizeFactors)
jg.controlResultsDeseq = results(jg.controlDeSeq)


# Run DeSeq on experiment
jg.experimentDeSeq<-jg.runDeSeq(jg.experimentPeaksetDeSeq, jg.conditions,jg.controlSizeFactors)
contrast <- c("Condition", "DOX", "none") ##specify base level for comparison
jg.experimentResultsDeseq   = results(jg.experimentDeSeq, contrast = contrast)
jg.experimentResultsDeseq_df <- as.data.frame(jg.experimentResultsDeseq)
#jg.experimentResultsDeseq_df_sig <- jg.experimentResultsDeseq_df #for saving raw DEseq2 output
jg.experimentResultsDeseq_df <- jg.experimentResultsDeseq_df[!is.na(jg.experimentResultsDeseq_df$padj),]
jg.experimentResultsDeseq_df_sig <- jg.experimentResultsDeseq_df[jg.experimentResultsDeseq_df$padj < 0.05 & 
                                                                   abs(jg.experimentResultsDeseq_df$log2FoldChange) > 0.5,]

##Peak annotation
sig_peak_bed <- do.call("rbind.data.frame", strsplit(gsub("-", ":", rownames(jg.experimentResultsDeseq_df_sig)), ":"))
colnames(sig_peak_bed) <- c("seqnames", "start", "end")
sig_peak_bed.GR = toGRanges(sig_peak_bed)
# Annotate peaks with information on closest TSS using precompiled annotation data
sig_peak_bed.anno =annotatePeakInBatch(sig_peak_bed.GR, AnnotationData=TSS.human.GRCh38)
##change settings to increase the span of peak location annotation as shown earlier line 93
library(org.Hs.eg.db)
sig_peak_bed.anno <- addGeneIDs(sig_peak_bed.anno, org.Hs.eg.db, mart = ensembl_hs, feature_id_type = "ensembl_gene_id",
                                IDs2Add = c("symbol", "entrezgene"))
sig_peak_bed.anno_df <- as.data.frame(sig_peak_bed.anno)

mapping_de <- select(org.Hs.eg.db,
                     keys = sig_peak_bed.anno_df$symbol,
                     columns = c("ENTREZID", "SYMBOL"),
                     keytype = "SYMBOL")
sig_peak_bed.anno_df$ENTREZID <- mapping[match(sig_peak_bed.anno_df$symbol, mapping$SYMBOL),2]

##Pathway analysis
library(ReactomePA)
pathways_peaks_de <- enrichPathway(sig_peak_bed.anno_df$ENTREZID)
dotplot(pathways_peaks_de)

library(enrichR)
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enrichedDE <- enrichr(sig_peak_bed.anno_df$symbol, dbs)
enrichedDE_MF <- enrichedDE[[1]]
enrichedDE_CC <- enrichedDE[[2]]
enrichedDE_BP <- enrichedDE[[3]]
write.table(enrichedDE_BP, "~/CHIP_seq/data_git/GO_BP_DEseq2_adj_pval_th05_GRCh38.tsv", sep = "\t", quote = F, row.names = F)


#Save results for example plots
#save(jg.experimentResultsDeseq,file="data/jg.experimentResultsDeseq.rda")
#save(jg.controlResultsDeseq,file="data/jg.controlResultsDeseq.rda")


## Example of jg.plotDeSeq
jg.plotDeSeq(jg.experimentResultsDeseq,
             p=0.01,
             title.main="Fold-change in TF binding",
             flip=T
)

##Combined figure

jg.plotDeSeqCombined(jg.controlResultsDeseq,
                     jg.experimentResultsDeseq,
                     title.main="TF Binding Folding changes on DOX treatment",
                     p=0.01, flip=TRUE)

