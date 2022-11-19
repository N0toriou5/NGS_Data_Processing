# TEPA Orazio
# concatenate together separate fastq files from the same sample:
cd /mnt/f/Projects/Orazio/TEPA/SAL9125/SAL9125A1
cat DIPG007_CON_1_S1_L001_R1_001.fastq.gz DIPG007_CON_1_S1_L002_R1_001.fastq.gz > /mnt/f/Projects/Orazio/TEPA/DIPG007_CON_1_S1_L001_R1_001.fastq.gz
cat DIPG007_CON_1_S1_L001_R2_001.fastq.gz DIPG007_CON_1_S1_L002_R2_001.fastq.gz > /mnt/f/Projects/Orazio/TEPA/DIPG007_CON_1_S1_L001_R2_001.fastq.gz
cd ..
rm -r SAL9125A1/

# fastqc
mkdir fastqc
for fname in *_001.fastq.gz
do
root=`basename $fname _001.fastq.gz`
echo "Doing $root"
fastqc ${root}_001.fastq.gz -o fastqc/ -t 6
done

# adapter trimming (not required in this case)
mkdir trimmed
for fname in *_R1_001.fastq.gz
do
root=`basename $fname _R1_001.fastq.gz`
echo "Doing $root"
fastp -i ${root}_R1_001.fastq.gz -I ${root}_R2_001.fastq.gz -o trimmed/${root}_R1_trimmed.fastq.gz -O trimmed/${root}_R2_trimmed.fastq.gz
done

# align with hisat2, point to hisat2 index
mkdir hisat_bams
index=/mnt/f/genomes/Human/hisat_index/hg38
# start alignment
for fname in *_R1_001.fastq.gz
do
cd /mnt/f/Projects/Orazio/TEPA/
root=`basename $fname _R1_001.fastq.gz`
echo "Doing $root"
hisat2 -x $index \
-p 6 \
-t \
-1 ${root}_R1_001.fastq.gz \
-2 ${root}_R2_001.fastq.gz \
| samtools view -bS - > hisat_bams/${root}.bam
# Sort and index BAMs
cd /mnt/f/Projects/Orazio/TEPA/hisat_bams
samtools sort -m 2G -@ 6 -O BAM -o ${root}.sorted.bam ${root}.bam 
samtools index ${root}.sorted.bam
rm ${root}.bam
done


### Convert to counts with the subread package
cd /mnt/f/Projects/Orazio/TEPA/hisat_bams
mkdir rawcounts
# humand gtf location: /mnt/f/genomes/Human/GTF
gtf=/mnt/f/genomes/Human/GTF/Homo_sapiens.GRCh38.103.gtf
featureCounts -T 7 -p -s 2 -t exon -g gene_name -a $gtf -o rawcounts/featureCounts.txt *sorted.bam
