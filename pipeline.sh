cd /bigdata/stajichlab/xxu092


git clone <http adress of our repository>
nano pipeline.sh


##download mrna.fasta of our species of interest, can use gzipped file for kallisto
curl -o mrna.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/324/715/GCF_003324715.1_ASM332471v1/GCF_003324715.1_ASM332471v1_rna_from_genomic.fna.gz

zmore mrna.fna.gz 

##download our SRA data using a loop 

nano download.sh
module load sratoolkit
##config https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration

#creat an accessionnumber list
nano sralist.tsv
#download raw reads
sbatch -a 1-12 download.sh

#run kallisto
sbatch -a 1-12 kallisto.sh

#write Rscript for DESeq2
nano kallisto_deseq.R

srun Rscript kallisto_deseq.R
