# the purpose of this script is to create a single count matrix for scATAC-seq data from species mapped to different genomes
# 
mkdir species_compare
cd species_compare
#
####################### 1) remove peaks that are unable to lift to either species genome ##############################
#
mkdir 0_liftPeaks
cd 0_liftPeaks
# 
##### a) species 1 (e.g. human)
#
mkdir human
cd human
#
# human consensus peaks (generated using the script scATAC_2_merge_singleCells.sh)
ln -s /directory_to_human_data/merge/consensus_peaks/consensus_peaks_summits_extend_nonOverlap.bed
#
## i) liftover from human -> chimp
liftover -minMatch=0.5 -bedPlus=7 consensus_peaks_summits_extend_nonOverlap.bed hg19ToPanTro4.over.chain peaks_hg19_to_panTro4.bed peaks_hg19_unmapped2panTro4.bed
#
# remove any peaks that lift to contigs/chrM/chrY 
bedtools intersect -a peaks_hg19_to_panTro4.bed -b /home/michael_boyle/1_Datasets/Genomes/panTro4/panTro4_chromSizes_noContigs.bed -wa | sort -k 1,1 -k 2,2n > peaks_hg19_to_panTro4_noContigs.bed
#
## ii) liftover from chimp -> human
#
liftover -minMatch=0.5 -bedPlus=7 peaks_hg19_to_panTro4_noContigs.bed panTro4ToHg19.over.chain peaks_hg19_to_panTro4_noContigs_back2species.bed peaks_hg19_to_panTro4_noContigs_back2species_unmapped.bed
#
# remove any peaks that lift to contigs/chrM/chrY 
bedtools intersect -a peaks_hg19_to_panTro4_noContigs_back2species.bed -b /home/michael_boyle/1_Datasets/Genomes/hg19/hg19_chromSizes_noContigs.bed -wa | sort -k 1,1 -k 2,2n > peaks_hg19_to_panTro4_noContigs_back2species_noContigs.bed
#
## iii) remove any peaks that don't overlap 50% of original peak 
bedtools intersect -a peaks_hg19_to_panTro4_noContigs_back2species_noContigs.bed -b consensus_peaks_summits_extend_nonOverlap.bed -wa -f 0.50 -r | sort -k 1,1 -k 2,2n -u > peaks_human_passFilter.bed
#
#
##### b) species 2 (e.g. chimp)
#
cd ..
mkdir chimp
cd chimp
#
# chimp consensus peaks (generated using the script scATAC_2_merge_singleCells.sh)
ln -s /directory_to_chimp_data/merge/consensus_peaks/consensus_peaks_summits_extend_nonOverlap.bed
#
## i) liftover from chimp -> human
liftover -minMatch=0.5 -bedPlus=7 consensus_peaks_summits_extend_nonOverlap.bed panTro4ToHg19.over.chain peaks_panTro4_to_hg19.bed peaks_panTro4_unmapped2hg19.bed
#
# remove any peaks that lift to contigs/chrM/chrY 
bedtools intersect -a peaks_panTro4_to_hg19.bed -b /home/michael_boyle/1_Datasets/Genomes/hg19/hg19_chromSizes_noContigs.bed -wa | sort -k 1,1 -k 2,2n > peaks_panTro4_to_hg19_noContigs.bed
#
## ii) liftover from human -> chimp
#
liftover -minMatch=0.5 -bedPlus=7 peaks_panTro4_to_hg19_noContigs.bed hg19ToPanTro4.over.chain peaks_panTro4_to_hg19_noContigs_back2species.bed peaks_panTro4_to_hg19_noContigs_back2species_unmapped.bed
#
# remove any peaks that lift to contigs/chrM/chrY 
bedtools intersect -a peaks_panTro4_to_hg19_noContigs_back2species.bed -b /home/michael_boyle/1_Datasets/Genomes/panTro4/panTro4_chromSizes_noContigs.bed -wa | sort -k 1,1 -k 2,2n > peaks_panTro4_to_hg19_noContigs_back2species_noContigs.bed
#
## iii) remove any peaks that don't overlap 50% of original peak 
bedtools intersect -a peaks_panTro4_to_hg19_noContigs_back2species_noContigs.bed -b consensus_peaks_summits_extend_nonOverlap.bed -wa -f 0.50 -r | sort -k 1,1 -k 2,2n -u > peaks_chimp_passFilter.bed
#
# iv) lift peaks that pass filter back to human genome so we can compare betweeen species (below)
liftover -minMatch=0.5 -bedPlus=7 peaks_chimp_passFilter.bed panTro4ToHg19.over.chain peaks_chimp_passFilter_hg19.bed peaks_chimp_passFilter_hg19_unmapped.bed
sort -k 1,1 -k 2,2n -u peaks_chimp_passFilter_hg19.bed | sponge peaks_chimp_passFilter_hg19.bed
#
#
####################### 2) create a merged peak list using peaks called in both species ##############################
#
cd ../../../species_compare
mkdir 1_mergePeaks
cd 1_mergePeaks
#
# symlink peaks that pass lift
ln -s /species_compare/human/peaks_human_passFilter.bed
ln -s /species_compare/chimp/peaks_chimp_passFilter_hg19.bed
#
# give each peak a name
awk 'OFS="\t"{print $1,$2,$3,"human_peak_"NR}' peaks_human_passFilter.bed | sort -k 1,1 -k 2,2n -u > human_passFilter_hg19.bed
awk 'OFS="\t"{print $1,$2,$3,"chimp_peak"NR}' peaks_chimp_passFilter_hg19.bed | sort -k 1,1 -k 2,2n -u > chimp_passFilter_hg19.bed
# 
# cat the peaks from each species
cat human_passFilter_hg19.bed chimp_passFilter_hg19.bed | sort -k 1,1 -k 2,2n > peaks_cat_hg19.bed
#
# merge the peaks
bedtools merge -i peaks_cat_hg19.bed -c 4 -o distinct | sort -k 1,1 -k 2,2n -u > peaks_cat_hg19_merge.bed
wc -l peaks_cat_hg19_merge.bed
#
# remove any peaks mapping to contigs/chrM/chrY (should be 0)
bedtools intersect -a peaks_cat_hg19_merge.bed -b /home/michael_boyle/1_Datasets/Genomes/hg19/hg19_chromSizes_noContigs.bed -wa | sort -k 1,1 -k 2,2n > peaks_cat_hg19_merge_noContigs.bed
wc -l peaks_cat_hg19_merge_noContigs.bed
#
# bring the merge file back to panTro4 coords (to create count matrix)
liftover -minMatch=0.5 -bedPlus=7 peaks_cat_hg19_merge_noContigs.bed hg19ToPanTro4.over.chain peaks_cat_hg19_noContigs_panTro4.bed peaks_cat_hg19_noContigs_panTro4_unmapped
wc -l peaks_cat_hg19_noContigs_panTro4.bed
#
# remove any peaks mapping to contigs/chrM/chrY (should be 0)
bedtools intersect -a peaks_cat_hg19_noContigs_panTro4.bed -b /home/michael_boyle/1_Datasets/Genomes/panTro4/panTro4_chromSizes_noContigs.bed -wa | sort -k 1,1 -k 2,2n > peaks_cat_hg19_noContigs_panTro4_noContigs.bed
#
#
## get peaks in SAF format for generating count matrix
awk 'OFS="\t"{print $4,$1,$2,$3,"-"}' peaks_cat_hg19_merge_noContigs.bed > peaks_cat_hg19_v1.SAF
awk 'OFS="\t"{print $4,$1,$2,$3,"-"}' peaks_cat_hg19_noContigs_panTro4_noContigs.bed > peaks_cat_panTro4_v1.SAF
#
# we need to rename the peaks to a single name, so we can merge the count matrices later on:
awk 'OFS="\t"{print $1,$2"_"$3"_"$4,"mergePeak_"NR}' peaks_cat_hg19_v1.SAF > peaks_cat_hg19_names.txt
#
# add the peak names
join -1 1 -2 1 -a1 -eNA -o 2.3,1.2,1.3,1.4,1.5 <(sort -k 1,1 peaks_cat_hg19_v1.SAF) <(sort -k 1,1 peaks_cat_hg19_names.txt) | sort -k 2,2 -k 3,3n | sed -e "s/\ /\t/g" > peaks_cat_hg19.SAF
join -1 1 -2 1 -a1 -eNA -o 2.3,1.2,1.3,1.4,1.5 <(sort -k 1,1 peaks_cat_panTro4_v1.SAF) <(sort -k 1,1 peaks_cat_hg19_names.txt) | sort -k 2,2 -k 3,3n | sed -e "s/\ /\t/g" > peaks_cat_panTro4.SAF
# 
awk 'OFS="\t"{print $2,$3,$1}' peaks_cat_hg19_names.txt | awk 'OFS="\t"{gsub("[_]","\t",$1)}1' > peaks_cat_hg19_names.bed
awk 'OFS="\t"{print $2,$3,$4,$1}' peaks_cat_panTro4.SAF > peaks_cat_panTro4_names.bed
#
#
####################### 3) generate count matrix using the merged peak list ##############################
#
cd ../../species_compare
mkdir 2_countMatrix
cd 2_countMatrix
#
#### a) human
# 
mkdir human
cd human
#
# create a link to the file that lists the directories to single cell BAM files and the header containing cell IDs (generated using the script scATAC_2_merge_singleCells.sh) 
ln -s /directory_to_human_data/merge/consensus_peaks/list_bamlist.txt 
ln -s /directory_to_human_data/merge/consensus_peaks/header.txt 
#
names=$(cat list_bamlist.txt)
#
featureCounts -p -a ../../1_mergePeaks/peaks_cat_hg19.SAF -F SAF -O -o matrix_human.txt $names
#
tail -n +3 matrix_human.txt | cut -f 1,7- > matrix_human_final.bed
#
cat header.txt matrix_human_final.bed > matrix_human_final_header.bed
#
#
#### b) chimp
#
cd ..
mkdir chimp
cd chimp
#
# create a link to the file that lists the directories to single cell BAM files and the header containing cell IDs (generated using the script scATAC_2_merge_singleCells.sh) 
ln -s /directory_to_chimp_data/merge/consensus_peaks/list_bamlist.txt 
ln -s /directory_to_chimp_data/merge/consensus_peaks/header.txt 
#
names=$(cat list_bamlist.txt)
#
featureCounts -p -a ../../1_mergePeaks/peaks_cat_panTro4.SAF -F SAF -O -o matrix_chimp.txt $names
#
tail -n +3 matrix_chimp.txt | cut -f 1,7- > matrix_chimp_final.bed
#
cat header.txt matrix_chimp_final.bed > matrix_chimp_final_header.bed

####################### 4) join the species in single count matrix and then binarize ##############################
#
cd ..
#
# do this in R (scATAC_analysis.R)
# see section 4) species comparisons (create 1 count matrix for all species)
#
###############################################################################################







