# Bonito
# https://github.com/nanoporetech/bonito/blob/master/README.md
# https://github.com/bioconvert/bioconvert/issues/88

Steps to Set Python3 as Default On ubuntu?
Check python version on terminal - python --version
Get root user privileges. On terminal type - sudo su
Write down the root user password.
Execute this command to switch to python 3.6. update-alternatives --install /usr/bin/python python /usr/bin/python3 1
Check python version - python --version
https://stackoverflow.com/questions/70588185/warning-the-script-pip3-8-is-installed-in-usr-local-bin-which-is-not-on-path
then install bonito with pip install ont-bonito

pip install -f https://download.pytorch.org/whl/torch_stable.html ont-bonito-cuda111
##
cd /mnt/f/Projects/IGA-Tech/data/pnigra
bonito basecaller dna_r10.4_e8.1_hac@v3.4 --recursive fast5/ > basecalls.fastq

PycoQC
### check nvidia toolkit version
sudo apt install nvidia-cuda-toolkit
nvcc --version

mkdir bonito
bonito basecaller dna_r10.4_e8.1_hac@v3.4 --recursive test/ > bonito/basecalls.fastq
#### reduced batchsize
bonito basecaller -v dna_r9.4.1_e8.1_hac@v3.3 --batchsize 12 --chunksize 1000 --overlap 500 --device cuda --recursive test/ > bonito/basecalls.fastq
bonito basecaller dna_r10.4_e8.1_hac@v3.4 --batchsize 12 --chunksize 10000 --overlap 500 --device cuda --recursive fast5/ > basecalls.fastq
### issue https://github.com/nanoporetech/bonito/issues/208

#### Cineca bonito

### install bonito developer into a venv 
https://github.com/nanoporetech/bonito
cd bonito
source venv3/bin/activate
bonito download --models --show #dna_r9.4.1_e8.1_hac@v3.3
#go to fast5 folder
cd /mnt/d/Projects/ONT/Data/pnigra/
bonito basecaller --device cuda dna_r9.4.1_e8.1_hac@v3.3 --recursive fast5/ > basecalls.fastq

### use not the dev version
conda activate
cd /mnt/d/Projects/ONT/Data/pnigra/
bonito basecaller --device cuda dna_r9.4.1_e8.1_hac@v3.3 --recursive fast5/ > basecalls.fastq

### create a reference file: https://simpsonlab-demo.readthedocs.io/en/latest/quickstart_eventalign.html
cd /mnt/d/Projects/ONT/Data/pnigra/

minimap2 -d ref.mmi ref.fa
reads=/mnt/d/Projects/ONT/Data/pnigra/fast5/
ref=/mnt/d/Projects/ONT/Data/pnigra/reference/pnigra.mmi
export OMP_NUM_THREADS=1 #https://github.com/nanoporetech/bonito/issues/140
bonito basecaller dna_r9.4.1_e8.1_hac@v3.3 --save-ctc --reference $ref $reads > /mnt/d/Projects/ONT/Data/pnigra/ctc-data/pnigra.sam
bonito train --directory /mnt/d/Projects/ONT/Data/pnigra/ctc-data/chunks.npy /data/training/model-dir


### Train a new Bonito model https://github.com/nanoporetech/bonito/issues/137
cd /mnt/d/Projects/ONT/Data/pnigra/reference/
minimap2 -x map-ont -d pnigra.mmi Pn_italica_V2_7112019.fasta
# index ground truth sequences
minimap2 -x map-ont -d ground_truth_reads.mmi ground_truth_reads.fa 
reads=/mnt/d/Projects/ONT/Data/pnigra/fast5/
ref=/mnt/d/Projects/ONT/Data/pnigra/reference/pnigra.mmi
dna=/home/notorious/bonito/bonito/data/dna_r9.4.1/
export OMP_NUM_THREADS=1
bonito basecaller dna_r9.4.1_e8_hac@v3.3 --reference $ref --save-ctc $reads > bonito_calls.sam  
