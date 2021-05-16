#!/bin/bash
# scATAC_2_merge_singleCells.sh
#_______________________________________________________________________________________________
#
# usage: 
# mkdir merge
# cd merge
# sh scATAC_2_merge_singleCells.sh bamlist.txt
#
#_______________________________________________________________________________________________
#
# The script first merges the processed single cells BAM files (generated with "scATAC_1_processing_singleCells.sh") into one merged BAM file using samtools
# The merged BAM is then used for peak calling with macs2
#_______________________________________________________________________________________________
#
# INPUT FILES REQUIRED
# $1: bamlist.txt = a list of directories of processed single cell BAM files (generated with the script: scATAC_1_processing_singleCells.sh) 
#_______________________________________________________________________________________________
#
# SOFTWARE TOOLS REQUIRED
# samtools - v1.3.1-21-g874baf3
# JAVA - Java-1.8
# PICARD - v2.1.1		# user needs to add directory to reference genome (e.g. /directory_to/Genomes/panTro4/whole_genome.fa)
# macs2 - v2.1.1.20160309
# bedtools - v2.25.0
# bedGraphToBigWig - v4
# liftOver   		      # user needs to add directory to the chain files in this script! (e.g. /directory_to/Genomes/liftOver/panTro4ToHg19.over.chain)
#_______________________________________________________________________________________________
#
# Written by: Michael James Boyle
# E-Mail: michael_boyle@eva.mpg.de
# Date Written: Aug-22-2017
###############################################################################################
#
# run the script
#
samtools merge scATAC_joint_noDup_noMT_MERGE.bam -b $1
samtools index scATAC_joint_noDup_noMT_MERGE.bam
#
# collect metrics
java -Xmx8g -jar picard.jar CollectMultipleMetrics I=scATAC_joint_noDup_noMT_MERGE.bam O=scATAC_joint_noDup_noMT_MERGE_metrics R=/directory_to/Genomes/panTro4/whole_genome.fa 
#
# call peaks
macs2 callpeak --nomodel --format BAMPE -t scATAC_joint_noDup_noMT_MERGE.bam -n scATAC_joint_noDup_noMT_MERGE --nolambda --keep-dup all --slocal 10000 --call-summits
#
cut -f 1,2,3,4,9,10 scATAC_joint_noDup_noMT_MERGE_peaks.narrowPeak | sort -k 1,1 -k 2,2n -k 5,5nr | sort -k 1,1 -k 2,2n -u > scATAC_joint_noDup_noMT_MERGE_peaks_narrowpeak_unique.bed
# 
# liftover if working with non-human primate date (ex. chimp -> human)
liftover -minMatch=0.5 scATAC_joint_noDup_noMT_MERGE_peaks_narrowpeak_unique.bed /directory_to/Genomes/liftOver/panTro4ToHg19.over.chain scATAC_joint_noDup_noMT_MERGE_peaks_narrowpeak_hg19.bed scATAC_joint_noDup_noMT_MERGE_peaks_narrowpeak_hg19_unmapped.bed
sort -k 1,1 -k 2,2n -u scATAC_joint_noDup_noMT_MERGE_peaks_narrowpeak_hg19.bed | sponge scATAC_joint_noDup_noMT_MERGE_peaks_narrowpeak_hg19.bed
#
###############
#summits
cut -f 1,2,3,4,9,10 scATAC_joint_noDup_noMT_MERGE_peaks.narrowPeak | sort -k 1,1 -k 2,2n | awk 'OFS="\t"{print $1,$2+$6,$2+$6+1,$4,$5}' - > scATAC_joint_noDup_noMT_MERGE_peaks_summits_MB.bed
awk 'OFS="\t"{print $1,$2-250,$2+250,$4,$5}' scATAC_joint_noDup_noMT_MERGE_peaks_summits_MB.bed | sort -k 1,1 -k 2,2n > scATAC_joint_noDup_noMT_MERGE_peaks_summits_extend.bed
# remove overlapping peaks by taking the more significant one
awk 'OFS="\t"{if (NR==1) {chr=$1; start=$2; end=$3; name=$4; score=$5} else {if ($1==chr && ($2>=start && $2<=end )) {if ($5>score) {chr=$1; start=$2; end=$3; name=$4; score=$5}} else {print chr,start,end,name,score; chr=$1; start=$2; end=$3; name=$4; score=$5} }}END{print chr,start,end,name,score}' scATAC_joint_noDup_noMT_MERGE_peaks_summits_extend.bed > scATAC_joint_noDup_noMT_MERGE_peaks_summits_extend_nonOverlap.bed
# get the top 50k peaks
sort -k 5,5nr scATAC_joint_noDup_noMT_MERGE_peaks_summits_extend_nonOverlap.bed | head -n 50000 > scATAC_joint_noDup_noMT_MERGE_peaks_summits_top50k.bed
sort -k 1,1 -k 2,2n scATAC_joint_noDup_noMT_MERGE_peaks_summits_top50k.bed > scATAC_joint_noDup_noMT_MERGE_peaks_accessibilityPeaks.bed
# liftOver summits
liftover -minMatch=0.5 scATAC_joint_noDup_noMT_MERGE_peaks_summits_extend_nonOverlap.bed /home/michael_boyle/1_Datasets/Genomes/liftOver/panTro4ToHg19.over.chain scATAC_joint_noDup_noMT_MERGE_peaks_summits_extend_nonOverlap_hg19.bed scATAC_joint_noDup_noMT_MERGE_peaks_summits_extend_nonOverlap_unmapped.bed 
##############
# Tn5 shift
samtools view -h scATAC_joint_noDup_noMT_MERGE.bam | awk 'BEGIN{FS=OFS="\t"}{if ($0 ~ /^@/) {print $0} else {if ($9>0) {print $1,$2,$3,$4+4,$5,$6,$7,$8,$9,$10,$11} else {if ($9<0) print $1,$2,$3,$4-5,$5,$6,$7,$8,$9,$10,$11}}}' | samtools view -b > scATAC_joint_noDup_noMT_MERGE_shift.bam
samtools sort -o scATAC_joint_noDup_noMT_MERGE_shift.sorted.bam scATAC_joint_noDup_noMT_MERGE_shift.bam
samtools index scATAC_joint_noDup_noMT_MERGE_shift.sorted.bam
rm scATAC_joint_noDup_noMT_MERGE_shift.bam
##############
# create bigWig
bedtools genomecov -ibam scATAC_joint_noDup_noMT_MERGE.bam -bg -g /directory_to/Genomes/panTro4/panTro4.chromSizes > scATAC_joint_noDup_noMT_MERGE.bg
sort -k 1,1 -k 2,2n scATAC_joint_noDup_noMT_MERGE.bg | sponge scATAC_joint_noDup_noMT_MERGE.bg
/directory_to/bedGraphToBigWig scATAC_joint_noDup_noMT_MERGE.bg /home/michael_boyle/1_Datasets/Genomes/panTro4/panTro4.chromSizes scATAC_joint_noDup_noMT_MERGE.bw 
#
#############
# get consensus peaks
mkdir consensus_peaks
cd consensus_peaks
sed ':a;N;$!ba;s/\n/ /g' ../$1 > list_bamlist.txt
names=$(cat list_bamlist.txt)
bedtools intersect -a ../scATAC_joint_noDup_noMT_MERGE_peaks_narrowpeak_unique.bed -b $names -wa -wb -filenames | cut -f 1,2,3,4,5,6,7 | bedtools merge -i stdin -c 4,5,6,7 -o distinct,distinct,distinct,count_distinct > cells_overlapping_peaks.bed
num_cells=$(cat ../$1 | wc -l)
num_cells2=$(($num_cells*5/100 ))
awk -v cell_min="$num_cells2" 'OFS="\t"{if ($7 >= cell_min) {print $1,$2,$3,$4,$5,$6,$7}}' cells_overlapping_peaks.bed | sort -k 1,1 -k 2,2n -u > consensus_peaks.bed
# intersect with summits
bedtools intersect -a ../scATAC_joint_noDup_noMT_MERGE_peaks_summits_extend_nonOverlap.bed -b consensus_peaks.bed -wa | sort -k 1,1 -k 2,2n -u > consensus_peaks_summits_extend_nonOverlap.bed
#############
# generate count matrix
#convert peaks to SAF format:
awk 'OFS="\t"{print $1"_"$2"_"$3,$1,$2,$3,"-"}' consensus_peaks_summits_extend_nonOverlap.bed > merge_consensus_peaks_summits_extend_nonOverlap.SAF
names=$(cat list_bamlist.txt)
featureCounts -p -a merge_consensus_peaks_summits_extend_nonOverlap.SAF -F SAF -O -o matrix.txt $names
tail -n +3 matrix.txt | cut -f 1,7- > matrix_final.bed
# create HEADER: 
sed -e 's/_scATAC_hg19_noDup_noMT.bam//g' -e 's/\/directory\/to\/singleCell\///g' ../$1 | sed -e 's/\//\t/g' | awk 'OFS="\t"{print $2}' | tr '\n' '\t' > header.txt
# add "peak" to first line of header.txt
sed -i 's/^/peak\t/g' header.txt
# add the header to matrix 
cat header.txt matrix_final.bed > matrix_final_header.bed
#
# create matrix for Cicero input
awk 'OFS="\t"{if (NR==1) {for(i=2;i<=NF;i++) header[i]=$i} else {for(i=2;i<=NF;i++){ print $1,header[i],$i}}}' matrix_final_header.bed > matrix_for_Cicero.bed
#
#############
done
exit 
#
