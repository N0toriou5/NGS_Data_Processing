# ABCC3 snippet code for TWIST1 ChIP-seq

cd /mnt/h/ABCC3
mkdir ChIP
# Download dataset from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80154

prefetch -c SRR3356388	#https://www.ncbi.nlm.nih.gov/sra/?term=SRR3356388
prefetch -c SRR3356385 	#https://www.ncbi.nlm.nih.gov/sra/?term=SRR3356385

# convert to fastq
fastq-dump -I --split-files SRR3356388 --outdir /mnt/h/ABCC3/ChIP # split each read in separate files (TWIST1)
fastq-dump -I --split-files SRR3356385 --outdir /mnt/h/ABCC3/ChIP # (Input)
fastqc SRR3356385_1.fastq
fastqc SRR3356388_1.fastq

# adapter trimming
mkdir trimmed/
fastp -i SRR3356385_1.fastq -o trimmed/input_trimmed.fastq
fastp -i SRR3356388_1.fastq -o trimmed/chip_trimmed.fastq
# FASTQC analysis
# high quality input fastq files, truseq adapters in it, decided not to trim the, (sw skipping)
# sequence duplication detected in twist1 files

# cutadapt -a AACCGGTT -o output.fastq input.fastq general command to use cutadapt

# using trim galore: 
#trim=/home/notorious/TrimGalore-0.6.0/trim_galore
#$trim --fastqc -j 5 SRR3356385_1.fastq
#$trim --fastqc -j 5 SRR3356388_1.fastq

#Input filename: SRR3356388_1.fastq
#Trimming mode: single-end
#Trim Galore version: 0.6.0
#Cutadapt version: 1.15
#Python version: 3.6.8
#Number of cores used for trimming: 5
#Quality Phred score cutoff: 20
#Quality encoding type selected: ASCII+33
#Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
#Maximum trimming error rate: 0.1 (default)
#Minimum required adapter overlap (stringency): 1 bp
#Minimum required sequence length before a sequence gets removed: 20 bp
#Running FastQC on the data once trimming has completed


# assign index variable
index=/mnt/f/genomes/Human/Index/Bowtie2/GRCh38_noalt_as/GRCh38_noalt_as
mkdir bowtie
# Align reads with bowtie2
bowtie2 -k 1 -p 6 -x $index -U input_trimmed.fastq | samtools view -uS - > bowtie/input.bam
bowtie2 -k 1 -p 6 -x $index -U chip_trimmed.fastq | samtools view -uS - > bowtie/twist.bam

# Sort and index BAMs
cd bowtie
samtools sort -m 2G -@ 6 -O BAM -o input.sorted.bam input.bam
samtools index input.sorted.bam
samtools sort -m 2G -@ 6 -O BAM -o twist.sorted.bam twist.bam
samtools index twist.sorted.bam

# Markup and remove duplicates

 PicardCommandLine MarkDuplicates \
 REMOVE_DUPLICATES=TRUE \
 I="twist.sorted.bam" O=twist_NODUPS.bam \
 M=twist_dup_metrics.txt 

 PicardCommandLine MarkDuplicates \
 REMOVE_DUPLICATES=TRUE \
 I="input.sorted.bam" O=input_NODUPS.bam \
 M=input_dup_metrics.txt 
 
# Blacklist regions

# ENCODE blacklists (save and unzip hg38)
# from http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/


# Blacklist with ENCODE regions
bl=/mnt/f/genomes/Human/blacklist/hg38.blacklist.bed

# Filter out black lists
mkdir filtered
bedtools intersect -v -abam input_NODUPS.bam -b $bl > filtered/input.filtered.bam
bedtools intersect -v -abam twist_NODUPS.bam -b $bl > filtered/twist.filtered.bam
samtools index filtered/input.filtered.bam
samtools index filtered/twist.filtered.bam


#Call MACS peaks (one to one samples)
cd filtered
macs2 callpeak -t twist.filtered.bam \
-c input.filtered.bam \
-f BAM -g hs -n twist -B -p 1e-9 --outdir macs_twist --verbose 2 --bdg 

# Prepare a bed file centered on transcript TSSs
gff3=/mnt/f/genomes/Human/annotations/gencode.v36.annotation.gff3
# Select only transcripts
grep -P "\ttranscript\t" $gff3 > alltranscripts_coords.bed
cut -f3 alltranscripts_coords.bed | sort | uniq

cd /mnt/h/ABCC3/ChIP
# Remove redundant genes in custom .bed file
uniq up.bed > tmp && mv tmp genesup.bed
less -S genesup.bed | wc -l
uniq dn.bed > tmp && mv tmp genesdn.bed
uniq nonhit.bed > tmp && mv tmp genesnot.bed
# Starting the heatmap plot on TSS (center) <- center on TSS signals

# Files needed:
### BAM from ChIP-Seq (TWIST1)
### BAM from input
### BED file with TSS coordinates

cd /mnt/h/ABCC3/ChIP/trimmed/bowtie/filtered
# Calculate number of aligned reads
samtools flagstat "input.filtered.bam"  # 26877009
samtools flagstat "twist.filtered.bam" # 15748498

# Convert BAM to BigWig
bamCoverage -p 6 -b "input.filtered.bam"  -o input.bw --scaleFactor 1
bamCoverage -p 6 -b "twist.filtered.bam" -o twist.bw --scaleFactor 1.70664 # That is, 26877009/15748498

# Subtract input from signal
bigwigCompare -b1 twist.bw -b2 input.bw --operation subtract -p 4 -o twist_diff.bw

# Choose reference point
bedup=/mnt/d/FABIT/Projects/ABCC3/genesup.bed
beddn=/mnt/d/FABIT/Projects/ABCC3/genesdn.bed
bednot=/mnt/d/FABIT/Projects/ABCC3/genesnot.bed
computeMatrix reference-point -p4 -S twist_diff.bw -R $bedup -a 10000 -b 10000 --referencePoint center -o upanalysis.mat.gz
computeMatrix reference-point -p4 -S twist_diff.bw -R $beddn -a 10000 -b 10000 --referencePoint center -o dnanalysis.mat.gz
computeMatrix reference-point -p4 -S twist_diff.bw -R $bednot -a 10000 -b 10000 --referencePoint center -o notanalysis.mat.gz


# Plot
plotHeatmap -m upanalysis.mat.gz -out plots/001_deeptools_up.png --averageTypeSummaryPlot mean \
--heatmapHeight 13 --heatmapWidth 6 --whatToShow "heatmap and colorbar" \
--colorList 'navajowhite,orange,black' \
--refPointLabel "TWIST1 peak centers" --plotTitle "ABCC3 targets vs. TWIST1 ChIP-Seq" \
--xAxisLabel "Distance from TSS center (bp)" --yAxisLabel "peaks" \
--dpi 600

plotHeatmap -m dnanalysis.mat.gz -out plots/001_deeptools_dn.png --averageTypeSummaryPlot mean \
--heatmapHeight 13 --heatmapWidth 6 --whatToShow "heatmap and colorbar" \
--colorList 'navajowhite,orange,black' \
--refPointLabel "TWIST1 peak centers" --plotTitle "ABCC3 targets vs. TWIST1 ChIP-Seq" \
--xAxisLabel "Distance from TSS center (bp)" --yAxisLabel "peaks" \
--dpi 600

plotHeatmap -m notanalysis.mat.gz -out plots/001_deeptools_not.png --averageTypeSummaryPlot mean \
--heatmapHeight 13 --heatmapWidth 6 --whatToShow "heatmap and colorbar" \
--colorList 'navajowhite,orange,black' \
--refPointLabel "TWIST1 peak centers" --plotTitle "ABCC3 targets vs. TWIST1 ChIP-Seq" \
--xAxisLabel "Distance from TSS center (bp)" --yAxisLabel "peaks" \
--dpi 600