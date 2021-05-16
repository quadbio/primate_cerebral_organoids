####################################################################################################################
#### script to analyze scATAC-seq data in R (using human data as example)
####################################################################################################################
#
# 1) chromVAR (motif, k-mer enrichment)
# 2) cicero (cluster cells) 
# 3) diffusion Map (obtain pseudotime) 
# 4) species comparisons (create 1 count matrix for all species)
#
####################################################################################################################
##################                1) chromVAR (motif, k-mer enrichment)
####################################################################################################################
#
setwd("/directory_to/chromVAR_analysis")
#
### input files needed:
# 1) peaks  (generated using the script "scATAC_2_merge_singleCells.sh")

# 2) file containing directories to processed BAM files (e.g. RG#_scATAC_hg19_noDup_noMT.bam) 
       # single cell BAM files were generated using the script "scATAC_1_processing_singleCells.sh"
       # the list of processed single cell BAMs is generated when using "scATAC_2_merge_singleCells.sh"

# 3) cell annotations (colData): each single cell has info with it's BAM file directory, an ID/name, cell line, age of the organoid, replicate number 
#
### load libraries
#
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(JASPAR2016)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
register(MulticoreParam(8))
BiocParallel::register(BiocParallel::SerialParam())
library(nabor)
library(BSgenome.Hsapiens.UCSC.hg19)
library(readr)
library(DBI)
data("human_pwms_v1") 
#
#
######### 1)  read in PEAKS, BAM files, cell annotations to generate count matrix ############
#
setwd("/directory_to/chromVAR_analysis/motif")
#
#### load PEAK file
peakfile <- "/directory_to/merge/consensus_peaks/cat_consensus_peaks_summits_extend_nonOverlap.bed"  # see script scATAC_2_merge_singleCells.sh
peaks <- getPeaks(peakfile, sort=TRUE)
#
### load BAM files 
bamfiles <- readLines("/directory_to/merge/bamlist.txt")  # for use with script scATAC_2_merge_singleCells.sh
#
## load cell annotations ("colData", contains )
library(readr)
colData <- read_delim("/directory_to/merge/colData.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
colData
#
### get the count matrix
#
fragment_counts <- getCounts(bamfiles, peaks, 
                              paired =  TRUE, 
                              by_rg = FALSE, 
                              format = "bam", 
                              colData)
#
####### 2)  add GC bias #######
#
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg19)
#
####### 3)  filter cells and peaks  #######
#
# filter cells based on min number of reads and nim # of reads in peaks (FRiP)
counts_filtered <- filterSamples(fragment_counts, min_depth = 5000,min_in_peaks = 0.05, shiny = F)
#
# plot the filters
filtering_plot <- filterSamplesPlot(fragment_counts, min_depth = 5000, min_in_peaks = 0.05, use_plotly = FALSE)
#
pdf(file = "filter_plot.pdf")
filtering_plot
dev.off()
#
#### filter peaks (remove overlapping)
counts_filtered <- filterPeaks(counts_filtered)
#
##################################################
### calculate deviations for each motif
# matchMotifs =  function matchMotifs from the motifmatchr package finds which peaks contain which motifs. 
motifs_2<-human_pwms_v1
motif_ix <- matchMotifs(motifs_2, counts_filtered,genome = BSgenome.Hsapiens.UCSC.hg19)
bg <- getBackgroundPeaks(object = counts_filtered)
dev <- computeDeviations(object = counts_filtered, 
                          annotations = motif_ix, 
                          background_peaks = bg)
#
####### Visualizing Deviations
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 30)
write.table(tsne_results,"tsne_results.txt",sep="\t",row.names=T,col.names = T)
tsne_results2<-tsne_results
colnames(tsne_results2)<-c("V2","V3")
#
###### Generate Plots 
#
tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "NEUROD6", 
                                 sample_column = "batch", shiny = F)

pdf(file = "cells_by_batch.pdf")
tsne_plots[[1]]
dev.off()
#
######
#
tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "NEUROD6", 
                                 sample_column = "depth", shiny = F)

pdf(file = "cells_by_depth.pdf")
tsne_plots[[1]]
dev.off()
#
######
# variability
# The first application is simply to compute the variability of each motif or annotation across the cells or samples of interest
variability <- computeVariability(dev)
#
pdf(file = "variability.pdf")
plotVariability(variability, use_plotly = F) 
dev.off()
#
sort_variability<-variability[order(variability$p_value),]
write.table(sort_variability,"variability.txt",sep="\t",row.names=T,col.names = T,quote=F)
#
#######################################
########### plot marker genes
#
##### SOX2
#
tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "SOX2", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "SOX2.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["SOX2"]][["data"]]
write.table(enrichment,"SOX2_enrichment.txt",sep="\t",row.names=T,col.names = T)

SOX2_enrichment<-enrichment["color"]
colnames(SOX2_enrichment)<-"SOX2"

annotation_motif<-cbind(tsne_results2,SOX2_enrichment)
plot<-ggplot(data=annotation_motif, aes(x=V2, y=V3, fill=SOX2))
plot<-plot+theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot<-plot+geom_point(shape=21, size=4)+scale_fill_gradient2(low="blue", high="red", breaks=c(-3,0,3), name="SOX2 Enrichment")+theme(legend.position="top")
ggsave(filename="SOX2_enrichment.pdf", width=5.6, height=4.9)

## SOX9

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "SOX9", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "SOX9.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["SOX9"]][["data"]]
write.table(enrichment,"SOX9_enrichment.txt",sep="\t",row.names=T,col.names = T)

SOX9_enrichment<-enrichment["color"]
colnames(SOX9_enrichment)<-"SOX9"

## EMX1

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "EMX1", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "EMX1.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["EMX1"]][["data"]]
write.table(enrichment,"EMX1_enrichment.txt",sep="\t",row.names=T,col.names = T)

EMX1_enrichment<-enrichment["color"]
colnames(EMX1_enrichment)<-"EMX1"



##########  NEUROD6

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "NEUROD6", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "NEUROD6.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["NEUROD6"]][["data"]]
write.table(enrichment,"NEUROD6_enrichment.txt",sep="\t",row.names=T,col.names = T)

NEUROD6_enrichment<-enrichment["color"]
colnames(NEUROD6_enrichment)<-"NEUROD6"

## NEUROD1
tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "NEUROD1", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "NEUROD1.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["NEUROD1"]][["data"]]
write.table(enrichment,"NEUROD1_enrichment.txt",sep="\t",row.names=T,col.names = T)

NEUROD1_enrichment<-enrichment["color"]
colnames(NEUROD1_enrichment)<-"NEUROD1"

## NHLH1
tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "NHLH1", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "NHLH1.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["NHLH1"]][["data"]]
write.table(enrichment,"NHLH1_enrichment.txt",sep="\t",row.names=T,col.names = T)

NHLH1_enrichment<-enrichment["color"]
colnames(NHLH1_enrichment)<-"NHLH1"

## EOMES
tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "EOMES", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "EOMES.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["EOMES"]][["data"]]
write.table(enrichment,"EOMES_enrichment.txt",sep="\t",row.names=T,col.names = T)

EOMES_enrichment<-enrichment["color"]
colnames(EOMES_enrichment)<-"EOMES"


################ get motif enrichment results

tsne_plots <- plotDeviationsTsne(dev, tsne_results, sample_column = "cell_state", shiny = F)

motif_zscores<-t(deviationScores(dev))
write.table(motif_zscores,"motif_zscores.txt",sep="\t",row.names=T,col.names = T,quote = F)

############### covariable motifs
# Covariability' is defined as covariance between Z-scores divided by variance of Z-scores for one motif

motif_covariable <- deviationsCovariability(dev)
write.table(motif_covariable,"motif_covariable.txt",sep="\t",row.names=T,col.names = T,quote = F)

############################## calculate  Differential accessibility 
# The differential_deviations function determines whether there is a significant difference between the bias corrected deviations for a given annotation between different groups. The groups can be specified by giving the column name of a column in the colData of the dev object or a vector of group assignments.
# see whether deviations differ between groups (deviations calculated above, which is used for tSNE)
# for example, we can see if the cells classified as "d60" and "d120" are signficiantly different from each other based on the deviation value
# there should BE difference in neurons vs progenitors, because they are enriched for different motifs
# we can add information about the assigned cell state to the colData cell annotations

diff_acc <- differentialDeviations(dev, "cell_state")
diff_acc<-diff_acc[order(diff_acc$p_value_adjusted),]
write.table(diff_acc,"diff_acc_cell_state.txt",sep="\t",row.names=T,col.names = T,quote=F)


############################## Differential variability
# The differential_variability function determines whether there is a significant difference between the variability of any of the annotations between different groups.
# determine whether groups differ in variability of a motif (for example, for NEUROD6)
# if we compare neurons to progenitors, then we hope to see a significant difference in the TFs that mark these different lines

diff_var <- differentialVariability(dev, "cell_state")
diff_var<-diff_var[order(diff_var$p_value_adjusted),]
write.table(diff_var,"diff_var_cell_state.txt",sep="\t",row.names=T,col.names = T,quote=F)


##########################################################################
#################### k-mers (7-mers)

setwd("/directory_to/chromVAR_analysis/kmer")

kmer_ix <- matchKmers(7, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)
kmer_dev <- computeDeviations(counts_filtered, kmer_ix, background_peaks = bg)

# plot k-mer results
tsne_results_kMer <- deviationsTsne(kmer_dev, threshold = 1.5, perplexity = 30)
write.table(tsne_results_kMer,"tsne_results_kMer.txt",sep="\t",row.names=T,col.names = T)
tsne_results_kMer2<-tsne_results_kMer
colnames(tsne_results_kMer2)<-c("V2","V3")


tsne_plots_kMer <- plotDeviationsTsne(kmer_dev, tsne_results_kMer, 
                                   sample_column = "day", shiny = F)

######

tsne_plots <- plotDeviationsTsne(kmer_dev, tsne_results_kMer,
                                 sample_column = "batch", shiny = F)

pdf(file = "kMer_batch.pdf")
tsne_plots[[1]]
dev.off()

######

tsne_plots <- plotDeviationsTsne(kmer_dev, tsne_results_kMer,
                                 sample_column = "depth", shiny = F)

pdf(file = "kMer_depth.pdf")
tsne_plots[[1]]
dev.off()


######

pdf(file = "kMer_CAGATGG_Ebox.pdf")
plotDeviationsTsne(kmer_dev, tsne_results_kMer, annotation = "CAGATGG", 
                     sample_column = "day", shiny = F)
dev.off()


###### plot motif enrichments on k-mer tSNE

# SOX9:
tsne_plots <- plotDeviationsTsne(dev, tsne_results_kMer, annotation = "SOX9", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "SOX9_kMer.pdf")
tsne_plots[[2]]
dev.off()

# NEUROD6
tsne_plots <- plotDeviationsTsne(dev, tsne_results_kMer, annotation = "NEUROD6", 
                                 sample_column = "day", shiny = FALSE)

pdf(file = "NEUROD6_kMer.pdf")
tsne_plots[[2]]
dev.off()


######
# calculate variability 
variability_kMer <- computeVariability(kmer_dev)

pdf(file = "variability_kMer.pdf")
plotVariability(variability_kMer, use_plotly = F) 
dev.off()

sort_variability_kmer<-variability_kMer[order(variability_kMer$p_value),]
write.table(sort_variability_kmer,"variability_kMer.txt",sep="\t",row.names=T,col.names = T,quote=F)


############### covariable kmers
# Covariability' is defined as covariance between Z-scores divided by variance of Z-scores for one motif

kmer_cov <- deviationsCovariability(kmer_dev)
write.table(kmer_cov,"kmer_covariable.txt",sep="\t",row.names=T,col.names = T,quote = F)

####### Differential accessibility (kMer)
# The differential_deviations function determines whether there is a significant difference between the bias corrected deviations for a given annotation between different groups. The groups can be specified by giving the column name of a column in the colData of the dev object or a vector of group assignments.

diff_acc_kmer <- differentialDeviations(kmer_dev, "cell_state")
diff_acc_kmer<-diff_acc_kmer[order(diff_acc_kmer$p_value_adjusted),]
write.table(diff_acc_kmer,"diff_acc_kmer_cell_state.txt",sep="\t",row.names=T,col.names = T,quote=F)

####### Differential variability (kMer)
# The differential_variability function determines whether there is a significant difference between the variability of any of the annotations between different groups.

diff_var_kmer <- differentialVariability(kmer_dev, "cell_state")
diff_var_kmer<-diff_var_kmer[order(diff_var_kmer$p_value_adjusted),]
write.table(diff_var_kmer,"diff_var_kmer_cell_state.txt",sep="\t",row.names=T,col.names = T,quote=F)


####### write k-mer scores to file
kmer_zscores<-t(deviationScores(kmer_dev))
write.table(kmer_zscores,"kmer_zscores.txt",sep="\t",row.names=T,col.names = T,quote = F)


###############################################################
#####################  other annotations beside motifs & k-mers

#### use our own annotations (ex. mouse scATAC-seq DA peaks per cluster, corresponding to different cell states)
### need to have all annotations in 1 bed file!

setwd("/directory_to/chromVAR_analysis/mouse_scATAC")

my_annotation_file <- '/directory_to/mouse_developing_brain/mouse_E11_P0_clusters_hg19.bed'

anno_ix <- getAnnotations(my_annotation_file, rowRanges = rowRanges(counts_filtered), column = 4)
head(anno_ix)

anno_colData <- read_delim("directory_to/mouse_developing_brain/anno_colData_names.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
anno_list<-as.list(anno_colData)
anno_ix@colData@listData<-anno_list
head(anno_ix)

######

dev_anno <- computeDeviations(object = counts_filtered, 
                              annotations = anno_ix, 
                              background_peaks = bg)

tsne_results_anno <- deviationsTsne(dev_anno, threshold = 1.5, perplexity = 30)

mousefetalBrain_zscores<-t(deviationScores(dev_anno))
write.table(mousefetalBrain_zscores,"mousefetalBrain_zscores.txt",sep="\t",row.names=T,col.names = T,quote = F)


####### Differential accessibility (mouse)

diff_acc_mouse <- differentialDeviations(dev_anno, "cell_state")
diff_acc_mouse<-diff_acc_mouse[order(diff_acc_mouse$p_value_adjusted),]
write.table(diff_acc_mouse,"diff_acc_mouse_cell_state.txt",sep="\t",row.names=T,col.names = T,quote=F)

####### Differential variability (mouse)
diff_var_mouse <- differentialVariability(dev_anno, "cell_state")
diff_var_mouse<-diff_var_mouse[order(diff_var_mouse$p_value_adjusted),]
write.table(diff_var_mouse,"diff_var_mouse_cell_state.txt",sep="\t",row.names=T,col.names = T,quote=F)


######################## plots (using tSNE coords from TFBS motif enrichment, above)

#####  Radial glia


tsne_plots <- plotDeviationsTsne(dev_anno, tsne_results, annotation = "K3", 
                                 sample_column = "batch", shiny = F)

pdf(file = "cells_by_K3_RG_motifs.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["K3"]][["data"]]
write.table(enrichment,"K3_enrichment.txt",sep="\t",row.names=T,col.names = T)


##### excitatory neurons

tsne_plots <- plotDeviationsTsne(dev_anno, tsne_results, annotation = "K13", 
                                 sample_column = "batch", shiny = F)

pdf(file = "cells_by_K13_excitatoryNeurons_motifs.pdf")
tsne_plots[[2]]
dev.off()

enrichment<-tsne_plots[["K13"]][["data"]]
write.table(enrichment,"K13_enrichment.txt",sep="\t",row.names=T,col.names = T)


######################## plots (using tSNE coords from 7-mer enrichment, above)

#####  Radial glia
tsne_plots <- plotDeviationsTsne(dev_anno, tsne_results_kMer, annotation = "K3", 
                                 sample_column = "batch", shiny = F)

pdf(file = "cells_by_K3_RG_kMer.pdf")
tsne_plots[[2]]
dev.off()


##### excitatory neurons

tsne_plots <- plotDeviationsTsne(dev_anno, tsne_results_kMer, annotation = "K13", 
                                 sample_column = "batch", shiny = F)

pdf(file = "cells_by_K13_excitatoryNeurons_kMer.pdf")
tsne_plots[[2]]
dev.off()


######################################
### make nice plots
library(ggplot2)

# KMER
setwd("/directory_to/chromVAR_analysis/nice_plots/kmer")
colnames(tsne_results_kMer)<-c("V2","V3")
annotation_row = colData(kmer_dev)[c("cell_state","day")]
colnames(annotation_row)<-c("cell_state","day")

annotation_rowEnd = colData(kmer_dev)[c("cell_state","cell","batch","depth")]
colnames(annotation_rowEnd)<-c("cell_name","cell","batch","depth")

TF_enrichments<-cbind(ARID5A_enrichment, ASCL1_enrichment, ASCL2_enrichment, ATF3_enrichment, ATOH1_enrichment, ATOH7_enrichment, ATOH8_enrichment, BACH1_enrichment, BACH2_enrichment, BATF_enrichment, BHLHA15_enrichment, BHLHE22_enrichment, BHLHE23_enrichment, DLX1_enrichment, DLX2_enrichment, DLX3_enrichment, DLX4_enrichment, DLX5_enrichment, DLX6_enrichment, ELK3_enrichment, EMX1_enrichment, EMX2_enrichment, EOMES_enrichment, FERD3L_enrichment, FIGLA_enrichment, FOS_enrichment, FOSB_enrichment, FOSL1_enrichment, FOXD4_enrichment, FOXJ1_enrichment, GBX1_enrichment, GBX2_enrichment, GSX1_enrichment, GSX2_enrichment, HESX1_enrichment, ID3_enrichment, ID4_enrichment, JUN_enrichment, JUNB_enrichment, JUND_enrichment, LEF1_enrichment, LHX1_enrichment, LHX2_enrichment, LHX3_enrichment, LHX4_enrichment, LHX5_enrichment, LHX6_enrichment, LHX9_enrichment, LMO2_enrichment, LYL1_enrichment, MAX_enrichment, MEOX1_enrichment, MESP1_enrichment, MYF5_enrichment, MYF5_enrichment, MYOD1_enrichment, NEUROD1_enrichment, NEUROD2_enrichment, NEUROD4_enrichment, NEUROD6_enrichment, NEUROG1_enrichment, NEUROG2_enrichment, NEUROG3_enrichment, NFIA_enrichment, NFIB_enrichment, NFIC_enrichment, NFIX_enrichment, NHLH1_enrichment, NHLH2_enrichment, NKX11_enrichment, NOTO_enrichment, NRF1_enrichment, OLIG1_enrichment, OLIG2_enrichment, OLIG3_enrichment, PAX4_enrichment, PAX6_enrichment, POU2F2_enrichment, POU2F3_enrichment, POU3F2_enrichment, POU3F4_enrichment, POU4F2_enrichment, POU5F1_enrichment, POU5F1B_enrichment, POU6F1_enrichment, POU6F2_enrichment, RFX1_enrichment, RFX2_enrichment, RFX3_enrichment, RFX5_enrichment, RFX6_enrichment, RFX7_enrichment, RFX8_enrichment, SCXA_enrichment, SCXB_enrichment, SHOX_enrichment, SMARCC1_enrichment, SNAI1_enrichment, SNAI2_enrichment, SOX10_enrichment, SOX2_enrichment, SOX5_enrichment, SOX6_enrichment, SOX9_enrichment, SRY_enrichment, TAL1_enrichment, TAL2_enrichment, TBR1_enrichment, TBX1_enrichment, TCF12_enrichment, TCF15_enrichment, TCF21_enrichment, TCF23_enrichment, TCF3_enrichment, TCF4_enrichment, TCF7L2_enrichment, TEAD1_enrichment, TEAD3_enrichment, TEAD4_enrichment, TFAP4_enrichment, TWIST2_enrichment, ZEB1_enrichment,PITX2_enrichment,ZBTB3_enrichment,GATA3_enrichment,NR2F2_enrichment,NR2F6_enrichment)
TF_short_enrichments<-cbind(ARID_enrichment, ASCL_enrichment, ATF_enrichment, ATOH_enrichment, BACH_enrichment, BATF_enrichment, BHLHE_enrichment, DLX_enrichment, ELK_enrichment, EMX_enrichment, FOS_enrichment, FOSB_enrichment, FOXD_enrichment, FOXJ_enrichment, GBX_enrichment, GSX_enrichment, HESX_enrichment, ID_enrichment, JUN_enrichment, JUNB_enrichment, JUND_enrichment, LEF_enrichment, LHX_enrichment, LMO_enrichment, LYL_enrichment, MAX_enrichment, MEOX_enrichment, MYF_enrichment, MYOD_enrichment, NEUROD_enrichment, NEUROG_enrichment, NFIA_enrichment_short, NFIB_enrichment_short, NFIC_enrichment_short, NFIX_enrichment_short, NHLH_enrichment, NKX_enrichment, NOTO_enrichment, NRF_enrichment, OLIG_enrichment, PAX_enrichment, POU_enrichment, RFX_enrichment, SHOX_enrichment, SOX_enrichment, SRY_enrichment, TAL_enrichment, TBR_enrichment, TCF_enrichment, TEAD_enrichment, TFAP_enrichment, TWIST_enrichment, ZEB_enrichment)
annos_all<-cbind(mousefetalBrain_zscores,adultBrain_zscores,DEgenes_zscores,DEgenes2_zscores,DEgenes3_zscores)

annotation_row2<-cbind(tsne_results_kMer,annotation_row,TF_enrichments,TF_short_enrichments,annos_all,annotation_rowEnd)
write.table(annotation_row2,"annos_all.txt",sep="\t",row.names=T,col.names = T,quote=F)

tsne_enrich_combine<-cbind(annotation_row2)
tsne_enrich_combine_df<-data.frame(tsne_enrich_combine)

plot<-ggplot(data=tsne_enrich_combine_df, aes(x=V2, y=V3, fill=as.factor(day)))
plot<-plot+theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot+geom_point(shape=21, size=4)+scale_fill_manual(values=c("red","orange","yellow","yellowgreen","blue","chartreuse3","purple"), name="Age")+theme(legend.position="top")
ggsave(filename="cells_by_day.pdf", width=5.6, height=4.9)

plot<-ggplot(data=tsne_enrich_combine_df, aes(x=V2, y=V3, fill=NEUROD6))
plot<-plot+theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank())
plot<-plot+geom_point(shape=21, size=4)+scale_fill_gradient2(low="blue", high="red", breaks=c(-3,0,3), name="NEUROD6 Enrichment")+theme(legend.position="top")
ggsave(filename="NEUROD6_enrichment.pdf", width=5.6, height=4.9)

#######################################

save.image("complete_datset.RData")


####################################################################################################################
###############                           2) cicero (cluster cells) 
####################################################################################################################
#
### input files needed:
# 1) count matrix for Cicero (matrix_for_Cicero.bed)     # generated with script "scATAC_2_merge_singleCells.sh"
#
# 2) cell annotations (colData): each single cell has info with it's BAM file directory, an ID/name, cell line, age of the organoid, replicate number 
#
# 3) peak annotations (peak_anno): metadata annotations for peaks, such as which ones fall in promoter/distal regions, nearest expressed TSS (within 1Mb)
#
### load libraries
#
library(cicero)
library(Sushi)
#
setwd("/directory_to/cicero_analysis")
#
#
###  import the data:
counts = read.table(file = "/directory_to/merge/consensus_peaks/matrix_for_Cicero.bed", header=F)  # generated with script "scATAC_2_merge_singleCells.sh"
input_cds <- make_atac_cds(counts, binarize = TRUE)
pData(input_cds) <- cbind(pData(input_cds), colData[row.names(pData(input_cds)),])
#
#
## add promoter info to file
peak_anno <- read.table(file = "/directory_to/peak_annos.tsv", header=T)
peak_anno <- data.frame(peak_anno)
input_cds <- annotate_cds_by_site(input_cds, peak_anno)
head(fData(input_cds))
colnames(fData(input_cds))[1]<-"site_name"
#
TSS_peak<-subset(fData(input_cds), !(canonical_promoter %in% NA))
TSS_peak<-row.names(TSS_peak)
#
##########################################################
# 1) aggregation of cell count in peaks 
#    parameters: distance between peaks
##########################################################
agg_cds <- aggregate_nearby_peaks(input_cds, distance = 0)  
agg_cds <- detectGenes(agg_cds)
agg_cds <- estimateSizeFactors(agg_cds)
agg_cds <- estimateDispersions(agg_cds)

# adding the reduced dimensions for each cell
redDim_coords<-tsne_results_kMer2
agg_cds@reducedDimA<-t(redDim_coords)
agg_cds <- clusterCells(agg_cds, verbose = F,num_clusters=2)

# plot clusters
png("plot_cell_clusters.png")
plot_cell_clusters(agg_cds, color_by = 'as.factor(Cluster)')
dev.off()

png("plot_cell_clusters_bySize_Factor.png")
plot_cell_clusters(agg_cds, color_by = 'Size_Factor')  + scale_colour_gradient2(low='blue', high = 'red')
dev.off()

png("plot_cell_clusters_by_num_genes_expressed.png")
plot_cell_clusters(agg_cds, color_by = 'num_genes_expressed') + scale_colour_gradient2(low='blue', high = 'red')
dev.off()

### identify DA peaks between clusters
clustering_DA_sites <- differentialGeneTest(agg_cds, #Takes a few minutes
                                            fullModelFormulaStr = '~Cluster', cores=32)

write.table(clustering_DA_sites,"clustering_DA_sites.txt",sep="\t",row.names=T,col.names = T)

# get the sites that are differentially accessible between clusters
clustering_DA_sites_signif <- rownames(clustering_DA_sites[which(clustering_DA_sites$qval<=0.05),])
write.table(clustering_DA_sites_signif,"clustering_DA_sites_signif_q0.05.txt",sep="\t",row.names=F,col.names = T,quote=F)

########## output cluster file names:
cluster1_cells<-pData(agg_cds[, pData(agg_cds)$Cluster %in% 1])$cells
write.table(cluster1_cells,"cluster1_cells.txt",sep="\t",row.names=F,col.names = F,quote=F)

cluster2_cells<-pData(agg_cds[, pData(agg_cds)$Cluster %in% 2])$cells
write.table(cluster2_cells,"cluster2_cells.txt",sep="\t",row.names=F,col.names = F,quote=F)

########## for each peak, calulcate the number of cells / total number of cells in cluster
# first, we need to get the count matrix back:
cds_exprs<-exprs(agg_cds)
cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
colnames(cds_exprs) <- c("f_id", "Cell", "expression")
cds_pData <- pData(agg_cds)[,c("Cluster","cell")]
cds_fData <- fData(agg_cds)
cds_exprs <- merge(cds_exprs, cds_fData, by.x="f_id", by.y="row.names")
cds_exprs <- merge(cds_exprs, cds_pData, by.x="Cell", by.y="row.names")

### now we need to calculate for each cluster 
library(data.table)

# cluster 1
cluster1_exprs<-subset(cds_exprs, Cluster %in% 1)
cluster1_exprs_df <- data.table(cluster1_exprs)
cluster1_exprs_df2 <- cluster1_exprs_df[,list(sumamount = sum(expression), normalize = ((sum(expression)/length(cluster1_cells))*100),normalize_allcells=((sum(expression)/num_cells_expressed)*100)), by = "f_id"]
cluster1_exprs_df3 <- cluster1_exprs_df2[,list(normalize_max = max(normalize)), by = "f_id"]
cluster1_exprs_df4 <- cluster1_exprs_df3[order(cluster1_exprs_df3$f_id)]
colnames(cluster1_exprs_df4)<-c("f_id","cluster1")

# cluster2
cluster2_exprs<-subset(cds_exprs, Cluster %in% 2)
cluster2_exprs_df <- data.table(cluster2_exprs)
cluster2_exprs_df2 <- cluster2_exprs_df[,list(sumamount = sum(expression), normalize = ((sum(expression)/length(cluster2_cells))*100),normalize_allcells=((sum(expression)/num_cells_expressed)*100)), by = "f_id"]
cluster2_exprs_df3 <- cluster2_exprs_df2[,list(normalize_max = max(normalize)), by = "f_id"]
cluster2_exprs_df4 <- cluster2_exprs_df3[order(cluster2_exprs_df3$f_id)]
colnames(cluster2_exprs_df4)<-c("f_id","cluster2")

merg1<-merge(cluster1_exprs_df4,cluster2_exprs_df4, by="f_id")
merg1[, "max"] <- apply(merg1[, 2:4], 1, max)
merg1[, "max_cluster"]<-colnames(merg1[,2:3])[apply(merg1[,2:3],1,which.max)]
write.table(merg1,"counts_across_clusters_allPeaks.txt",sep="\t",row.names=F,col.names = T,quote=F)
peak_assigned_cluster <- merg1[,c(1,6)]

####
# get significanlly D.A. peaks that overlap promoters 
signif_sites<-rownames(clustering_DA_sites[which(clustering_DA_sites$qval<=0.05),])
TSS_signif <-subset(signif_sites,signif_sites %in% TSS_peak)
TSS_signif_name<-subset(TSS, TSS$site_name %in% TSS_signif)
TSS_signif_name2<-merge(clustering_DA_sites,TSS_signif_name,by.x=0,by.y=0)
# add the cluster to it
TSS_signif_name3<-merge(TSS_signif_name2,peak_assigned_cluster,by.x="Row.names",by.y="f_id")
length(TSS_signif_name3)
write.table(TSS_signif_name3,"signif_sites_promoters.txt",sep="\t",row.names=F,col.names = T,quote=F)


#### lets get the top 250 peaks per cluster

signif_clusteringDAsites<-clustering_DA_sites[which(clustering_DA_sites$qval<=0.05),]
signif_clusteringDAsites2<-merge(signif_clusteringDAsites,peak_assigned_cluster,by.x="site_name",by.y="f_id")
rownames(signif_clusteringDAsites2)<-signif_clusteringDAsites2$site_name
write.table(signif_clusteringDAsites2,"signif_clusteringDAsites2.txt",sep="\t",row.names=F,col.names = T,quote=F)

# top 250:
cluster1_DApeaks<-subset(signif_clusteringDAsites2, max_cluster %in% "cluster1")
cluster1_DApeaks_250<-row.names(cluster1_DApeaks)[order(cluster1_DApeaks$qval)][1:250]
write.table(cluster1_DApeaks_250,"cluster1_DApeaks_250.txt",sep="\t",row.names=F,col.names = T,quote=F)

cluster2_DApeaks<-subset(signif_clusteringDAsites2, max_cluster %in% "cluster2")
cluster2_DApeaks_250<-row.names(cluster2_DApeaks)[order(cluster2_DApeaks$qval)][1:250]
write.table(cluster2_DApeaks_250,"cluster2_DApeaks_250.txt",sep="\t",row.names=F,col.names = T,quote=F)

top250_each_cluster<-c(cluster1_DApeaks_250,cluster2_DApeaks_250)
write.table(top250_each_cluster,"top250_each_cluster.txt",sep="\t",row.names=F,col.names = F,quote=F)

######## plot some of the peaks (cluster)
peaks_to_plot <-row.names(clustering_DA_sites)[order(clustering_DA_sites$qval)][1:10]
# specific_peaks_to_plot <- row.names(subset(fData(agg_cds),site_name %in% c("chr1_713871_714370", "chr1_762697_763196")))
plot_genes_jitter(agg_cds[peaks_to_plot,],grouping = "Cluster",color_by = "Cluster")

pdf("jitter_peaks_in_clusters.pdf") 
plot_genes_jitter(agg_cds[peaks_to_plot,],grouping = "Cluster",color_by = "Cluster")
dev.off()

## plot maker genes
TSS_marker2<-subset(TSS, canonical_promoter %in% c("NEUROD6", "EOMES", "SOX2", "TBR1", "PAX6"))
TSS_marker_peak2<-as.character(row.names(TSS_marker2))

pdf("jitter_peaks_in_clusters_markerGene_TSS.pdf") 
plot_genes_jitter(agg_cds[TSS_marker_peak2,],grouping = "Cluster",color_by = "Cluster")
dev.off()

########
final_results<-pData(agg_cds)
write.table(final_results,"final_results_aggCDS.txt",sep="\t",row.names=F,col.names = T,quote=F)

########
save.image("complete_datset.RData")


####################################################################################################################
###############                           3) diffusion Map (obtain pseudotime) 
####################################################################################################################

setwd("/directory_to/diffusion_analysis")

### input files needed:
# count matrix: cells as columns, peaks as rows 

### load libraries
#
library('Matrix')
library('parallel')
library('MASS')
library('diffusionMap')
library('FNN')
library('igraph')
library('princurve')
library('ggplot2')
library('inline')
library('gplots')
library('scatterplot3d')
#
# read in count matrx
counts<- read.delim("/directory_to/merge/counsensus_peaks/matrix_final_header.bed", header=T)    # generated with script "scATAC_2_merge_singleCells.sh"
rownames(counts)<-counts$peak 

########
counts_DA<-subset(counts, peak %in% row.names(top250_each_cluster))
nrow(counts_DA)
write.table(counts_DA,"counts_DA.tsv",sep="\t",row.names=F,col.names = T,quote=F)
#
# we need only the site_name and counts for the analysis
counts_DA_short<-counts_DA[c(1,8:ncol(counts_DA))]
# rownames(counts_DA_short)<-counts_DA_short$site_name
counts_DA_short<-counts_DA_short[2:ncol(counts_DA_short)]

# convert to matrix
matrx <-data.matrix(counts_DA_short,rownames.force = T)

# binarize matrix
bin_matrx <- ifelse(matrx>1,1,matrx)

# distance matrix: 1-Pearson_correlation (cell-cell)
dist_matrx <- 1 - cor(bin_matrx)

max_dim=50
eigen_values.red.th=0.04

# create diffusion map
diffusion_map <- diffuse(dist_matrx, neigen=max_dim, maxdim=max_dim)

# get eigenvalues for each Diffusion Component
eigen_values <- diffusion_map$eigenvals

# eigenvlaue / sum(eigenvalues)
eigen_values.red <- eigen_values/sum(eigen_values)

# take only those which are greater than 0.04 (4% of total values)
evdim <- rev(which(eigen_values.red > eigen_values.red.th))[1]
evdim <- max(2, evdim, na.rm=TRUE)
evdim
# if using 3 dimensions:
# evdim<-3

# get the diffusion coordinates for each cell
colnames(diffusion_map$X) <- paste0('DMC', 1:ncol(diffusion_map$X))
res <- diffusion_map$X[, 1:evdim]

# scale the coordinates 
res <- scale(diffusion_map$X[, 1:evdim])
rownames(res)<-rownames(dist_matrx)
write.table(res,"diffusion_map_coords.txt",sep="\t",row.names=T,col.names = T)

# organoid stage pseudotime: order along DC1
pseudotime_DC1<-as.data.frame(res[,"DMC1"])
colnames(pseudotime_DC1)<-"DC1"
pseudotime_DC1$rank <- NA
order.scores<-order(pseudotime_DC1$DC1,decreasing=T)
pseudotime_DC1$rank[order.scores] <- 1:nrow(pseudotime_DC1)
pseudotime_DC1<-pseudotime_DC1[order(pseudotime_DC1$DC1),]
pseudotime_DC1$rank<-pseudotime_DC1$rank/nrow(pseudotime_DC1)
pseudotime_DC1$cells<-rownames(pseudotime_DC1)
write.table(res,"pseudotime_DC1.txt",sep="\t",row.names=T,col.names = T)

# for whole trajectory data: repeat the above analysis to obtain diffusion map coordinates
# use principle curve to obtain pseudotime values and combine with the pseudotime values obtained for organoid stage

################# plot 3-D diffusion map (when using data for whole trajectory)
library(scatterplot3d)

plot3d <- scale(diffusion_map$X[, 1:3])
plot3d<-as.data.frame(plot3d)
plot3d_coords_names<-data.frame(plot3d,rownames(dist_matrx))
colnames(plot3d_coords_names)<-c("DMC1","DMC2","DMC3","cell_name")

plot3d_combine <- merge(plot3d_coords_names,colData,by.x="cell_name",by.y="cell", sort=F)

plot3d_combine$cell_state <- as.factor(plot3d_combine$cell_state)
scatterplot3d(plot3d_combine$DMC1,plot3d_combine$DMC2,plot3d_combine$DMC3, pch=16, highlight.3d=F, main="3D Scatterplot") 

##### add colors to cell state
plot3d_combine$cell_state <- factor(plot3d_combine$cell_state, levels = c("iPSC","EB","NEcto","NEpith","NPC","neuron"))

levels(plot3d_combine$cell_state)
colors <- c("#313695","#4575b4","deepskyblue2","aquamarine2","gold1","darkred")
colors_cell_state <- colors[as.numeric(plot3d_combine$cell_state)]


# plot in 3D
pdf("diffusion_plot_cellColors.pdf")
scatterplot3d(plot3d_combine$DMC1,plot3d_combine$DMC2,plot3d_combine$DMC3, pch=16, highlight.3d=F, main="Diffusion Map",xlab="DC1", ylab="DC2", zlab="DC3",color=colors_cell_state) 
legend(x=4,y=6, legend = levels(plot3d_combine$cell_state),col = colors, pch = 16, inset = -0.15, xpd = TRUE, horiz = F)
dev.off()

# plot in 2D:
# levels(plot3d_combine$cell_state)
plot<-ggplot(data=plot3d_combine, aes(x=DMC1, y=DMC2, fill=as.factor(cell_state)))
plot+geom_point(shape=21, size=4)+scale_fill_manual(values=c("red","orange","yellow","green","blue","violet"), name="cell state")+theme(legend.position="top")
ggsave(filename="Diffuse_by_cell_state.pdf", width=5.6, height=4.9)


#### draw principal curve through diffusion map  

# guide line by ordering from starting point (POU5F1 enriched cells)
pricu.f=1/3
ordercells <-res[with(res, order(-DMC1)),]
ordercells2<-data.matrix(ordercells)
pricu2 <- principal.curve(res[,1:3], smoother='lowess', start=ordercells2[1:20,],trace=TRUE, f=pricu.f, stretch=333)
pc.line2 <- as.data.frame(pricu2$s[order(pricu2$lambda), ])

pdf("princurve_DC_ordered.pdf")
plot(res)
lines(pricu2$s[ order (pricu2$lambda), ], lty = 1, lwd = 4, col = "purple", type = "l")
dev.off()


##### add princcurve to 3d plot:

pdf("3Dplot_with_trajectory_DC1_DC2_DC3_withLIne.pdf")
s3d<-scatterplot3d(plot3d_combine$DMC1,plot3d_combine$DMC2,plot3d_combine$DMC3, pch=21, highlight.3d=F, main="Diffusion Map",
                   xlab="DC1", ylab="DC2", zlab="DC3",color="black",box=FALSE,bg=colors_cell_state, cex.symbols = 1.75) 
s3d$points3d(pc.line2$DMC1,pc.line2$DMC2,pc.line2$DMC3,pch=24,type="l",lwd = 5)
dev.off()


#####  pseudotime value is assigned to each cell as projection along principal curve
pseudotime_princurve = pricu2$lambda
cells_assigned_pseudotime<-data.frame(res,pseudotime_princurve)
cells_assigned_pseudotime_names<- data.frame(cells_assigned_pseudotime,rownames(dist_matrx))
colnames(cells_assigned_pseudotime_names)<-c("DMC1","DMC2","DMC3","pseudotime","cell_name")
write.table(cells_assigned_pseudotime_names,"pseudotime_earlyTimePoints.txt",sep="\t",row.names=F,col.names = T,quote = F)
# add in pseudotime values obtained from diffusion map on organoid stage only, repeat plots above


########
save.image("complete_datset.RData")



####################################################################################################################
###############                           4) species comparisons (create 1 count matrix for all species)
####################################################################################################################
#
setwd("/directory_to/species_compare/2_countMatrix")  
#
### human count matrix:
counts_human = read.table(file = "matrix_human_final_header.bed", header=T)     # generated using "scATAC_3_merge_species.txt"
row.names(counts_human)<-counts_human$peak
nrow(counts_human)
ncol(counts_human)
#
### chimp count matrix:
counts_chimp = read.table(file = "matrix_chimp_final_header.bed", header=T)      # generated using "scATAC_3_merge_species.txt"
row.names(counts_chimp)<-counts_chimp$peak
nrow(counts_chimp)
ncol(counts_chimp)
#
####### merge human and chimp count matrices
merge_counts<-merge(counts_human,counts_chimp,by="peak")
nrow(merge_counts)
ncol(merge_counts)
row.names(merge_counts)<-merge_counts$peak
#
# rename the peaks by their shared name (mergePeak1,etc)
peak_names = read.table(file = "../1_mergePeaks/peaks_cat_hg19_names.txt", header=F)     # generated using "scATAC_3_merge_species.txt"
colnames(peak_names) <- c("peak","site_name","mergePeak_name")
nrow(peak_names)
# 
merge_counts_named<-merge(peak_names,merge_counts,by.x="mergePeak_name",by.y="peak")
nrow(merge_counts_named)
ncol(merge_counts_named)
#
write.table(merge_counts_named,"merge_counts_nonBinary.tsv",sep="\t",row.names=F,col.names = T,quote=F)
#
######## binarize the merged count matrix
# 
mat<-data.matrix(merge_counts,rownames.force = T)
# binarize
bmat <- ifelse(mat>1,1,mat)
merge_counts2<-as.data.frame(bmat)
merge_counts2$peak<-rownames(merge_counts2)
nrow(merge_counts2)
ncol(merge_counts2)
#
### add the peak names
merge_counts_named_binarize<-merge(peak_names,merge_counts2,by.x="mergePeak_name",by.y="peak")
nrow(merge_counts_named_binarize)
ncol(merge_counts_named_binarize)
#
############ 
# add peak annos
peak_annos = read.table(file = "../1_mergePeaks/peaks_cat_hg19_names.bed", header=T)
peak_annos<-peak_annos[c(1:4)]
#
merge_counts_annos<-merge(peak_annos,merge_counts_named_binarize,by.x="mergePeak_name",by.y="mergePeak_name")
#
merge_counts_annos<-merge_counts_annos[c(2:4,1,5:ncol(merge_counts_annos))]
colnames(merge_counts_annos)[1]<-"chr"
colnames(merge_counts_annos)[2]<-"start"
colnames(merge_counts_annos)[3]<-"end"
#
head(merge_counts_annos,1)
# chr start end peak_name cell1 cell2 cell3... 
#
# sort by chr then start
merge_counts_annos_sort<-merge_counts_annos[with(merge_counts_annos, order(chr, start)),]
#
########## write the final binary count matrix in BED format and as tsv
#
write.table(merge_counts_annos_sort,"merge_counts_annos.bed",sep="\t",row.names=F,col.names = T,quote=F)
#
write.table(merge_counts_annos_sort[c(4:ncol(merge_counts_annos_sort))],"merge_counts_annos.tsv",sep="\t",row.names=F,col.names = T,quote=F)
#
save.image("complete_datset.RData")
#
################################################################################################################