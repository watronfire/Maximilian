#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=1 -l mem=2gb -q workq -o /gpfs/home/natem/logs/Maximilian.txt -j oe

export BLASTDB=/gpfs/group/andersen/natem/db/blast_db

rm ~/logs/Maximilian.txt
cd /gpfs/home/natem/scripts/Maximilian
snakemake -k -j 50 --configfile=snakemake_config.json --cluster-config cluster.json --cluster "qsub -V -l walltime={cluster.walltime} -l mem={cluster.mem} -l nodes={cluster.n} -q {cluster.queue} -o {cluster.logfile} -j {cluster.stdout}"
