##Brundle workflow for sharing
##sample sheets Experiment_sample_sheet.csv and Control_sample_sheet.csv are shared on github data folder

library(Brundle)

##Initial conditions
jg.controlMinOverlap      <- 5
jg.treatedCondition       =  "DOX"
jg.untreatedCondition     =  "none"

jg.experimentSampleSheet = "/projects/drc_scratch/chipseq_wf/Chipseq_mycn/data/Experiment_sample_sheet.csv"
jg.controlSampleSheet    = "/projects/drc_scratch/chipseq_wf/Chipseq_mycn/data/Control_sample_sheet.csv"

##Create Diffbind object
dbaExperiment <- jg.getDba(jg.experimentSampleSheet, bRemoveDuplicates=TRUE)
dbaControl    <- jg.getDba(jg.controlSampleSheet, bRemoveDuplicates=TRUE)

saveRDS(dbaExperiment, "diffbind_exp", compress = T) ##Please share
saveRDS(dbaControl, "diffbind_cont", compress = T) ##Please share

##Read counts from peak data

jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
jg.controlPeakset    <- jg.dbaGetPeakset(dbaControl)

saveRDS(jg.experimentPeakset, "control_counts_expt_object", compress = T) ##Please share
saveRDS(jg.controlPeakset, "control_counts_control_object", compress = T) ##Please share


#Get counts for the treated control samples.
jg.controlCountsTreated<-jg.getControlCounts(jg.controlPeakset,
                                             jg.controlSampleSheet,
                                             jg.treatedCondition )
#Repeat for the untreated/control samples
jg.controlCountsUntreated<-jg.getControlCounts(jg.controlPeakset,
                                               jg.controlSampleSheet,
                                               jg.untreatedCondition)


saveRDS(jg.controlCountsTreated, "control_counts_treated_object",compress = T) ##Please share
saveRDS(jg.controlCountsUntreated, "control_counts_untreated_object", compress = T) ##Please share

##Calculation of correction factor
#Get the sample names for replicates that represent the two conditions.
jg.untreatedNames <- names(jg.controlCountsUntreated)
jg.treatedNames   <- names(jg.controlCountsTreated)

#Calculate correction factor for DiffBind
jg.correctionFactor<-jg.getCorrectionFactor(jg.experimentSampleSheet,
                                            jg.treatedNames,
                                            jg.untreatedNames)

saveRDS(jg.correctionFactor, "correctionFactor_object", compress = T) ##Please share
