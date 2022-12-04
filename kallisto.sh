#!/usr/bin/bash -l
#SBATCH -p short --mem 24gb -N 1 -n 16  --out logs/kallisto.%a.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

module load kallisto
kallisto index -i Rsphae.idx  mrna.fna.gz
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
     echo "cannot run without a number provided either cmdline or --array in sbatch"
     exit
 fi
fi

IFS=,
cat samplerun.csv | sed -n ${N}p | while read ACC GENO COND REP
do
 OUT=output/$GENO.$COND.$REP
 kallisto quant -t $CPU --single -l 300 -s 20 -i Rsphae.idx -o $OUT sra_data/${ACC}_1.fastq.gz
done
