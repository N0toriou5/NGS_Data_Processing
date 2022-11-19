### Test Bonito vs. Guppy with Vitis vinifera reads
phase 1: Guppy test
# https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html
FLO-MIN112 r10.4 pores
conda activate
cd /mnt/d/Projects/ONT/
guppy_basecaller --print_workflows
fast5=/mnt/d/Projects/ONT/Data/R10_vitis_vinifera/20220330_1450_MN28002_FAT11181_c81c15f5/fast5
guppy_basecaller -c dna_r10.4_e8.1_hac.cfg -i $fast5 -s guppy_basecall -x 'auto' --recursive --compress_fastq 
# Init time: 1249 ms
# Caller time: 4983466 ms, Samples called: 25058397260, samples/s: 5.02831e+06

conda deactivate
# quality check with pycoQC
~/.local/bin/pycoQC -f /mnt/d/Projects/ONT/Data/R10_vitis_vinifera/20220330_1450_MN28002_FAT11181_c81c15f5/guppy_basecall/sequencing_summary.txt --html_outfile local_run.html

### cat all passed fastqs
cd pass
cat *.fastq.gz > all_guppy.fastq.gz


# read mapping
conda activate
index=/mnt/d/Projects/ONT/Data/R10_vitis_vinifera/reference/Vitis_vinifera.PN40024.v4.dna.toplevel.fa.gz
fastq=/mnt/d/Projects/ONT/Data/R10_vitis_vinifera/20220330_1450_MN28002_FAT11181_c81c15f5/guppy_basecall/pass/all_guppy.fastq.gz
mkdir aligned
minimap2 $index $fastq > aligned/Vvinifera.sam -ax map-ont
cd aligned
samtools view -bS pnigra.sam > pnigra.bam
samtools sort pnigra.bam -o pnigra.sorted.bam
samtools index pnigra.sorted.bam

## check alignment metrics
picard CollectAlignmentSummaryMetrics \
          R=$index \
          I=pnigra.sorted.bam \
          O=metrics.txt



