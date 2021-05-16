#!/bin/bash
# scATAC_1_processing_singleCells.sh
#_______________________________________________________________________________________________
# usage: 
# cd /directory/with_single_cell_fastq_folders
# sh scATAC_1_processing_singleCells.sh /directory/with_single_cell_fastq_folders *
#
#_______________________________________________________________________________________________
#
# The script first aligns NGS reads against an indexed reference assembly using bowtie2 v2.2.9 
# The alignment BAMs produced above are processed using samtools and the PICARD pipeline
# The processed BAMs are then used for peak calling using macs2
#_______________________________________________________________________________________________
#
# INPUT FILES REQUIRED
# $1 --> single cell folder directory, where each folder represents a single cell (defined by it's read group ID ["RG#")] and contains fastq files) (e.g. human_day60_RG1)
# $i: RG# (not REQURIED as input by user; input via bin/ls)
#_______________________________________________________________________________________________
#
# SOFTWARE TOOLS REQUIRED
# bowtie2 - v2.2.9 	      # user needs to put path to the genome in the script (e.g. /directory_to/Genomes/hg19/bowtie2/hg19)
# JAVA - Java-1.8
# samtools - v1.3.1-21-g874baf3
# PICARD - v2.1.1		# user needs to add directory to reference genome (e.g. /directory_to/Genomes/hg19/whole_genome.fa)
# macs2 - v2.1.1.20160309 
# bedtools - v2.25.0
# bedGraphToBigWig - v4
#_______________________________________________________________________________________________
#
# Written by: Michael James Boyle
# E-Mail: michael_boyle@eva.mpg.de
# Date Written: Aug-22-2017
###############################################################################################
#
# run the script
#
for i in `/bin/ls -d *RG*` ; do
   	cd $1/$i
#
# align to hg19 (human reference genome)
	nohup bowtie2 -X 2000 -p 12 -x /directory_to/Genomes/hg19/bowtie2/hg19 -1 "$i"_r1.fastq -2 "$i"_r2.fastq -S "$i"_scATAC_hg19.sam
	mkdir hg19
	mv *nohup* hg19/
	# keep paired reads, sort, and index
	samtools view -f 2 -b -S -o "$i"_scATAC_hg19.bam "$i"_scATAC_hg19.sam
	samtools sort "$i"_scATAC_hg19.bam -o "$i"_scATAC_hg19_sorted.bam
	samtools index "$i"_scATAC_hg19_sorted.bam
	# mark duplicates (Picard tools)
	java -Xmx8g -jar picard.jar MarkDuplicates I="$i"_scATAC_hg19_sorted.bam OUTPUT="$i"_scATAC_hg19_noDup.bam ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
	samtools index "$i"_scATAC_hg19_noDup.bam
	# retain reads with MAPQ30 ; remove chrM and chrY
	samtools view -b -h -q 30 "$i"_scATAC_hg19_noDup.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX > "$i"_scATAC_hg19_noDup_q30.bam
	samtools sort "$i"_scATAC_hg19_noDup_q30.bam -o "$i"_scATAC_hg19_noDup_q30_sorted.bam
	samtools index "$i"_scATAC_hg19_noDup_q30_sorted.bam
	# remove reads overlapping the blacklist (chrM nuclear homologs)
	samtools view "$i"_scATAC_hg19_noDup_q30_sorted.bam -b -h -L /directory_to/mitochondrial_blacklist.hg19_sort.bed -o "$i"_scATAC_hg19_noDup_MT.bam -U "$i"_scATAC_hg19_noDup_noMT.bam
	samtools index "$i"_scATAC_hg19_noDup_noMT.bam
	# collect metrics
	java -Xmx8g -jar picard.jar CollectMultipleMetrics I="$i"_scATAC_hg19_noDup_noMT.bam O="$i"_scATAC_hg19_metrics_noMT R=/directory_to/Genomes/hg19/whole_genome.fa 
	java -Xmx8g -jar picard.jar CollectAlignmentSummaryMetrics I="$i"_scATAC_hg19_sorted.bam O="$i"_scATAC_hg19_alignment_metrics_sorted.txt R=/directory_to/Genomes/hg19/whole_genome.fa
	# estimate library complexity 
	java -Xmx8g -jar picard.jar EstimateLibraryComplexity I="$i"_scATAC_hg19_noDup_noMT.bam O="$i"_scATAC_hg19_libComplex_metrics.txt
	# call peaks 
	macs2 callpeak --nomodel --format BAMPE -t "$i"_scATAC_hg19_noDup_noMT.bam -n "$i"_scATAC_hg19_noDup_noMT --nolambda --keep-dup all --slocal 10000 --call-summits
	cut -f 1,2,3 "$i"_scATAC_hg19_noDup_noMT_peaks.narrowPeak | uniq > "$i"_scATAC_hg19_noDup_noMT_narrowPeak_GREAT.bed
	cut -f 1,2,3,4,9,10 "$i"_scATAC_hg19_noDup_noMT_peaks.narrowPeak | sort -k 1,1 -k 2,2n -k 5,5nr | sort -k 1,1 -k 2,2n -u > "$i"_scATAC_hg19_noDup_noMT_narrowPeak_peaks.bed
	# get summits
	awk 'OFS="\t"{print $1,$2+$6,$2+$6+1,$4,$5}' "$i"_scATAC_hg19_noDup_noMT_narrowPeak_peaks.bed | sort -k 1,1 -k 2,2n > "$i"_scATAC_hg19_summits.bed
	# extended summits
	awk 'OFS="\t"{print $1,$2-250,$2+249,$4,$5}' "$i"_scATAC_hg19_summits.bed | sort -k 1,1 -k 2,2n > "$i"_scATAC_hg19_summits_extend.bed
	# remove overlapping peaks
	awk 'OFS="\t"{if (NR==1) {chr=$1; start=$2; end=$3; name=$4; score=$5} else {if ($1==chr && ($2>=start && $2+25<=end )) {if ($5>score) {chr=$1; start=$2; end=$3; name=$4; score=$5}} else {print chr,start,end,name,score; chr=$1; start=$2; end=$3; name=$4; score=$5} }}END{print chr,start,end,name,score}' "$i"_scATAC_hg19_summits_extend.bed > "$i"_scATAC_hg19_summits_extend_nonOverlap.bed
	# take the top 50,000
	sort -k 5,5nr "$i"_scATAC_hg19_summits_extend_nonOverlap.bed | head -n 50000 > "$i"_scATAC_hg19_summits_top50k.bed
	sort -k 1,1 -k 2,2n "$i"_scATAC_hg19_summits_top50k.bed > "$i"_scATAC_hg19_summits_top50k_accessibilityPeaks.bed
	# get matrix for Cicero
	awk -v cell="$i" 'OFS="\t"{print $1"_"$2"_"$3,cell,"1"}' "$i"_scATAC_hg19_summits_extend.bed > "$i"_scATAC_hg19_summits_extend_forCicero.bed
	# create Tn5 shift file
	samtools view -h "$i"_scATAC_hg19_noDup_noMT.bam | awk 'BEGIN{FS=OFS="\t"}{if ($0 ~ /^@/) {print $0} else {if ($9>0) {print $1,$2,$3,$4+4,$5,$6,$7,$8,$9,$10,$11} else {if ($9<0) print $1,$2,$3,$4-5,$5,$6,$7,$8,$9,$10,$11}}}' | samtools view -b > "$i"_scATAC_hg19_noDup_noMT_shift.bam
	samtools sort -o "$i"_scATAC_hg19_noDup_noMT_shift.sorted.bam "$i"_scATAC_hg19_noDup_noMT_shift.bam
	samtools index "$i"_scATAC_hg19_noDup_noMT_shift.sorted.bam
	rm "$i"_scATAC_hg19_noDup_noMT_shift.bam
	# create bigWig file
	bedtools genomecov -ibam "$i"_scATAC_hg19_noDup_noMT.bam -bg -g /directory_to/Genomes/hg19/chrom.sizes > "$i"_scATAC_hg19_noDup_noMT.bg
	sort -k 1,1 -k 2,2n "$i"_scATAC_hg19_noDup_noMT.bg | sponge "$i"_scATAC_hg19_noDup_noMT.bg
	bedGraphToBigWig "$i"_scATAC_hg19_noDup_noMT.bg /directory_to/Genomes/hg19/chrom.sizes "$i"_scATAC_hg19_noDup_noMT.bw
	#clean up
	rm "$i"_scATAC_hg19.sam
	rm "$i"_scATAC_hg19_summits.bed
	rm "$i"_scATAC_hg19_summits_extend.bed
	mv *bam hg19/
	mv *bai hg19/
	mv *bed hg19/
	mkdir mets
	mv *metrics* mets/
	mv mets/ hg19/
	mv *xls hg19/
	mv *peaks* hg19/
	mkdir wiggle
	mv *bg wiggle/
	mv *bw wiggle/
	mv wiggle/ hg19/
#
# now for pantro4 (chimpanzee reference genome)
#
	nohup bowtie2 -X 2000 -p 12 -x /directory_to/Genomes/panTro4/bowtie2/panTro4 -1 "$i"_r1.fastq -2 "$i"_r2.fastq -S "$i"_scATAC_panTro4.sam
	mkdir panTro4
	mv *nohup* panTro4/
	# keep paired reads, sort, and index
	samtools view -f 2 -b -S -o "$i"_scATAC_panTro4.bam "$i"_scATAC_panTro4.sam
	samtools sort "$i"_scATAC_panTro4.bam -o "$i"_scATAC_panTro4_sorted.bam
	samtools index "$i"_scATAC_panTro4_sorted.bam
	# mark duplicates (Picard tools)
	java -Xmx8g -jar picard.jar MarkDuplicates I="$i"_scATAC_panTro4_sorted.bam OUTPUT="$i"_scATAC_panTro4_noDup.bam ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1023
	samtools index "$i"_scATAC_panTro4_noDup.bam
	# retain reads with MAPQ30 ; remove chrM and chrY
	samtools view -b -h -q 30 "$i"_scATAC_panTro4_noDup.bam chr1 chr2A chr2B chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX > "$i"_scATAC_panTro4_noDup_q30.bam
	samtools sort "$i"_scATAC_panTro4_noDup_q30.bam -o "$i"_scATAC_panTro4_noDup_q30_sorted.bam
	samtools index "$i"_scATAC_panTro4_noDup_q30_sorted.bam
	# remove reads overlapping the blacklist (chrM nuclear homologs)
	samtools view "$i"_scATAC_panTro4_noDup_q30_sorted.bam -b -h -L /directory_to/mitochondrial_blacklist.panTro4_sort.bed -o "$i"_scATAC_panTro4_noDup_MT.bam -U "$i"_scATAC_panTro4_noDup_noMT.bam
	samtools index "$i"_scATAC_panTro4_noDup_noMT.bam
	# collect metrics
	java -Xmx8g -jar picard.jar CollectMultipleMetrics I="$i"_scATAC_panTro4_noDup_noMT.bam O="$i"_scATAC_panTro4_metrics_noMT R=/directory_to/Genomes/panTro4/whole_genome.fa
	java -Xmx8g -jar picard.jar CollectAlignmentSummaryMetrics I="$i"_scATAC_panTro4_sorted.bam O="$i"_scATAC_panTro4_metrics_sortedBam.txt R=/directory_to/Genomes/panTro4/whole_genome.fa
	# estimate library complexity 
	java -Xmx8g -jar picard.jar EstimateLibraryComplexity I="$i"_scATAC_panTro4_noDup_noMT.bam O="$i"_scATAC_panTro4_noDup_noMT_libComplex_metrics.txt
	# call peaks
	macs2 callpeak --nomodel --format BAMPE -t "$i"_scATAC_panTro4_noDup_noMT.bam -n "$i"_scATAC_panTro4_noDup_noMT --nolambda --keep-dup all --slocal 10000 --call-summits
	cut -f 1,2,3 "$i"_scATAC_panTro4_noDup_noMT_peaks.narrowPeak | uniq > "$i"_scATAC_panTro4_noDup_noMT_narrowPeak.bed
	cut -f 1,2,3,4,9,10 "$i"_scATAC_panTro4_noDup_noMT_peaks.narrowPeak | sort -k 1,1 -k 2,2n -k 5,5nr | sort -k 1,1 -k 2,2n -u > "$i"_scATAC_panTro4_peaks_narrowpeak.bed
	# liftOver
	liftover -minMatch=0.5 "$i"_scATAC_panTro4_peaks_narrowpeak.bed /directory_to/Genomes/liftOver/panTro4ToHg19.over.chain "$i"_scATAC_panTro4_peaks_narrowpeak_hg19.bed "$i"_scATAC_panTro4_unMapped_peaks 
	# get summits
	awk 'OFS="\t"{print $1,$2+$6,$2+$6+1,$4,$5}' "$i"_scATAC_panTro4_peaks_narrowpeak.bed | sort -k 1,1 -k 2,2n > "$i"_scATAC_panTro4_summits.bed
	# extend summits
	awk 'OFS="\t"{print $1,$2-250,$2+249,$4,$5}' "$i"_scATAC_panTro4_summits.bed | sort -k 1,1 -k 2,2n > "$i"_scATAC_panTro4_summits_extend.bed
	# remove overlapping peaks
	awk 'OFS="\t"{if (NR==1) {chr=$1; start=$2; end=$3; name=$4; score=$5} else {if ($1==chr && ($2>=start && $2+25<=end )) {if ($5>score) {chr=$1; start=$2; end=$3; name=$4; score=$5}} else {print chr,start,end,name,score; chr=$1; start=$2; end=$3; name=$4; score=$5} }}END{print chr,start,end,name,score}' "$i"_scATAC_panTro4_summits_extend.bed > "$i"_scATAC_panTro4_summits_extend_nonOverlap.bed
	# liftover 
	liftover -minMatch=0.5 "$i"_scATAC_panTro4_summits_extend_nonOverlap.bed /directory_to/Genomes/liftOver/panTro4ToHg19.over.chain "$i"_scATAC_panTro4_summits_extend_nonOverlap_hg19.bed "$i"_scATAC_panTro4_nonOverlap_unMapped_peaks 
	# take the top 50,000
	sort -k 5,5nr "$i"_scATAC_panTro4_summits_extend_nonOverlap.bed | head -n 50000 > "$i"_scATAC_panTro4_summits_top50k.bed
	sort -k 1,1 -k 2,2n "$i"_scATAC_panTro4_summits_top50k.bed > "$i"_scATAC_panTro4_accessibilityPeaks.bed
	# get matrix for Cicero
	awk -v cell="$i" 'OFS="\t"{print $1"_"$2"_"$3,cell,"1"}' "$i"_scATAC_panTro4_summits_extend.bed > "$i"_scATAC_panTro4_summits_extend_forCicero.bed
	# create Tn5 shift file
	samtools view -h "$i"_scATAC_panTro4_noDup_noMT.bam | awk 'BEGIN{FS=OFS="\t"}{if ($0 ~ /^@/) {print $0} else {if ($9>0) {print $1,$2,$3,$4+4,$5,$6,$7,$8,$9,$10,$11} else {if ($9<0) print $1,$2,$3,$4-5,$5,$6,$7,$8,$9,$10,$11}}}' | samtools view -b > "$i"_scATAC_panTro4_noDup_noMT_shift.bam
	samtools sort -o "$i"_scATAC_panTro4_noDup_noMT_shift.sorted.bam "$i"_scATAC_panTro4_noDup_noMT_shift.bam
	samtools index "$i"_scATAC_panTro4_noDup_noMT_shift.sorted.bam
	rm "$i"_scATAC_panTro4_noDup_noMT_shift.bam
	# create bigWig
	bedtools genomecov -ibam "$i"_scATAC_panTro4_noDup_noMT.bam -bg -g /directory_to/Genomes/panTro4/panTro4.chromSizes > "$i"_scATAC_panTro4_noDup_noMT.bg
	sort -k 1,1 -k 2,2n "$i"_scATAC_panTro4_noDup_noMT.bg | sponge "$i"_scATAC_panTro4_noDup_noMT.bg
	bedGraphToBigWig "$i"_scATAC_panTro4_noDup_noMT.bg /directory_to/Genomes/panTro4/panTro4.chromSizes "$i"_scATAC_panTro4_noDup_noMT.bw 
	# clean up
	rm "$i"_scATAC_panTro4.sam
	rm "$i"_scATAC_panTro4_summits.bed
	rm "$i"_scATAC_panTro4_summits_extend.bed
	mv *bam panTro4/
	mv *bai panTro4/
	mv *bed panTro4/
	mkdir mets
	mv *metrics* mets/
	mv mets/ panTro4/
	mv *xls panTro4/
	mv *peaks* panTro4/
	mkdir wiggle
	mv *bg wiggle/
	mv *bw wiggle/
	mv wiggle/ panTro4/
#
# now for rheMac8 (macaque reference genome)
#
	nohup bowtie2 -X 2000 -p 12 -x /directory_to/Genomes/rheMac8/bowtie2/rheMac8 -1 "$i"_r1.fastq -2 "$i"_r2.fastq -S "$i"_scATAC_rheMac8.sam
	mkdir rheMac8
	mv *nohup* rheMac8/
	# keep paired reads, sort, and index
	samtools view -f 2 -b -S -o "$i"_scATAC_rheMac8.bam "$i"_scATAC_rheMac8.sam
	samtools sort "$i"_scATAC_rheMac8.bam -o "$i"_scATAC_rheMac8_sorted.bam
	samtools index "$i"_scATAC_rheMac8_sorted.bam
	# mark duplicates (Picard tools)
	java -Xmx8g -jar picard.jar MarkDuplicates I="$i"_scATAC_rheMac8_sorted.bam OUTPUT="$i"_scATAC_rheMac8_noDup.bam ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1023
	samtools index "$i"_scATAC_rheMac8_noDup.bam
	# retain reads with MAPQ30 ; remove chrM and chrY
	samtools view -b -h -q 30 "$i"_scATAC_rheMac8_noDup.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX > "$i"_scATAC_rheMac8_noDup_q30.bam
	samtools sort "$i"_scATAC_rheMac8_noDup_q30.bam -o "$i"_scATAC_rheMac8_noDup_q30_sorted.bam
	samtools index "$i"_scATAC_rheMac8_noDup_q30_sorted.bam
	# remove reads overlapping the blacklist (chrM nuclear homologs)
	samtools view "$i"_scATAC_rheMac8_noDup_q30_sorted.bam -b -h -L /directory_to/mitochondrial_blacklist.rheMac8_sort.bed -o "$i"_scATAC_rheMac8_noDup_MT.bam -U "$i"_scATAC_rheMac8_noDup_noMT.bam
	samtools index "$i"_scATAC_rheMac8_noDup_noMT.bam
	# collect metrics
	java -Xmx8g -jar picard.jar CollectMultipleMetrics I="$i"_scATAC_rheMac8_noDup_noMT.bam O="$i"_scATAC_rheMac8_metrics_noMT R=/directory_to/Genomes/rheMac8/rheMac8.fa
	java -Xmx8g -jar picard.jar CollectAlignmentSummaryMetrics I="$i"_scATAC_rheMac8_sorted.bam O="$i"_scATAC_rheMac8_metrics_sortedBam.txt R=/directory_to/Genomes/rheMac8/rheMac8.fa
	# estimate library complexity 
	java -Xmx8g -jar picard.jar EstimateLibraryComplexity I="$i"_scATAC_rheMac8_noDup_noMT.bam O="$i"_scATAC_rheMac8_noDup_noMT_libComplex_metrics.txt
	# call peaks
	macs2 callpeak --nomodel --format BAMPE -t "$i"_scATAC_rheMac8_noDup_noMT.bam -n "$i"_scATAC_rheMac8_noDup_noMT --nolambda --keep-dup all --slocal 10000 --call-summits
	cut -f 1,2,3 "$i"_scATAC_rheMac8_noDup_noMT_peaks.narrowPeak | uniq > "$i"_scATAC_rheMac8_noDup_noMT_narrowPeak.bed
	cut -f 1,2,3,4,9,10 "$i"_scATAC_rheMac8_noDup_noMT_peaks.narrowPeak | sort -k 1,1 -k 2,2n -k 5,5nr | sort -k 1,1 -k 2,2n -u > "$i"_scATAC_rheMac8_peaks_narrowpeak.bed
	# liftOver
	liftover -minMatch=0.5 "$i"_scATAC_rheMac8_peaks_narrowpeak.bed /directory_to/Genomes/liftOver/rheMac8ToHg19.over.chain "$i"_scATAC_rheMac8_peaks_narrowpeak_hg19.bed "$i"_scATAC_rheMac8_unMapped_peaks 
	# get summits
	awk 'OFS="\t"{print $1,$2+$6,$2+$6+1,$4,$5}' "$i"_scATAC_rheMac8_peaks_narrowpeak.bed | sort -k 1,1 -k 2,2n > "$i"_scATAC_rheMac8_summits.bed
	# extend summits
	awk 'OFS="\t"{print $1,$2-250,$2+249,$4,$5}' "$i"_scATAC_rheMac8_summits.bed | sort -k 1,1 -k 2,2n > "$i"_scATAC_rheMac8_summits_extend.bed
	# remove overlapping peaks
	awk 'OFS="\t"{if (NR==1) {chr=$1; start=$2; end=$3; name=$4; score=$5} else {if ($1==chr && ($2>=start && $2+25<=end )) {if ($5>score) {chr=$1; start=$2; end=$3; name=$4; score=$5}} else {print chr,start,end,name,score; chr=$1; start=$2; end=$3; name=$4; score=$5} }}END{print chr,start,end,name,score}' "$i"_scATAC_rheMac8_summits_extend.bed > "$i"_scATAC_rheMac8_summits_extend_nonOverlap.bed
	# liftover
	liftover -minMatch=0.5 "$i"_scATAC_rheMac8_summits_extend_nonOverlap.bed /directory_to/Genomes/liftOver/rheMac8ToHg19.over.chain "$i"_scATAC_rheMac8_summits_extend_nonOverlap_hg19.bed "$i"_scATAC_rheMac8_nonOverlap_unMapped_peaks 
	# take the top 50,000
	sort -k 5,5nr "$i"_scATAC_rheMac8_summits_extend_nonOverlap.bed | head -n 50000 > "$i"_scATAC_rheMac8_summits_top50k.bed
	sort -k 1,1 -k 2,2n "$i"_scATAC_rheMac8_summits_top50k.bed > "$i"_scATAC_rheMac8_accessibilityPeaks.bed
	# get matrix for Cicero
	awk -v cell="$i" 'OFS="\t"{print $1"_"$2"_"$3,cell,"1"}' "$i"_scATAC_rheMac8_summits_extend.bed > "$i"_scATAC_rheMac8_summits_extend_forCicero.bed
	# create Tn5 shift file
	samtools view -h "$i"_scATAC_rheMac8_noDup_noMT.bam | awk 'BEGIN{FS=OFS="\t"}{if ($0 ~ /^@/) {print $0} else {if ($9>0) {print $1,$2,$3,$4+4,$5,$6,$7,$8,$9,$10,$11} else {if ($9<0) print $1,$2,$3,$4-5,$5,$6,$7,$8,$9,$10,$11}}}' | samtools view -b > "$i"_scATAC_rheMac8_noDup_noMT_shift.bam
	samtools sort -o "$i"_scATAC_rheMac8_noDup_noMT_shift.sorted.bam "$i"_scATAC_rheMac8_noDup_noMT_shift.bam
	samtools index "$i"_scATAC_rheMac8_noDup_noMT_shift.sorted.bam
	rm "$i"_scATAC_rheMac8_noDup_noMT_shift.bam
	# create bigWig
	bedtools genomecov -ibam "$i"_scATAC_rheMac8_noDup_noMT.bam -bg -g /directory_to/Genomes/rheMac8/rheMac8.chromSizes > "$i"_scATAC_rheMac8_noDup_noMT.bg
	sort -k 1,1 -k 2,2n "$i"_scATAC_rheMac8_noDup_noMT.bg | sponge "$i"_scATAC_rheMac8_noDup_noMT.bg
	bedGraphToBigWig "$i"_scATAC_rheMac8_noDup_noMT.bg /directory_to/Genomes/rheMac8/rheMac8.chromSizes "$i"_scATAC_rheMac8_noDup_noMT.bw 
	# clean up
	rm "$i"_scATAC_rheMac8.sam
	rm "$i"_scATAC_rheMac8_summits.bed
	rm "$i"_scATAC_rheMac8_summits_extend.bed
	mv *bam rheMac8/
	mv *bai rheMac8/
	mv *bed rheMac8/
	mkdir mets
	mv *metrics* mets/
	mv mets/ rheMac8/
	mv *xls rheMac8/
	mv *peaks* rheMac8/
	mkdir wiggle
	mv *bg wiggle/
	mv *bw wiggle/
	mv wiggle/ rheMac8/
#
#
done
