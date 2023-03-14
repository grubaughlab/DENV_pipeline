#!/usr/bin/env bash

module load dSQ

cp /gpfs/ycga/project/grubaugh/shared/DENVSEQ/SCRIPT/REF_PRIMERS/* .;

rm -rf jobs.txt >/dev/null 2>&1; 

cat $1 | while read i; 
    do 
        echo "bash /gpfs/ycga/project/grubaugh/shared/DENVSEQ/SCRIPT/DENV_MAPPER.sh ${i}/*/${i}*_R1_*.fastq.gz ${i}/*/${i}*_R2_*.fastq.gz DENV.refs.txt" >> jobs.txt; 
    done

rm -rf DENV.serotype.calls.tsv >/dev/null 2>&1

dsqCmd=`echo -e $(dsq --job-name denv.mapper --job-file jobs.txt --mem-per-cpu=10G --cpus-per-task=1) | tr ":" "\n" | grep sbatch`;
dsqJobID=$(eval $dsqCmd | awk '{print $4}');

sbatch --job-name=denv.summary --depend=afterok:${dsqJobID} /gpfs/ycga/project/grubaugh/shared/DENVSEQ/SCRIPT/DENV_summarise.sh >/dev/null 2>&1;

rm -rf dsq-jobs-* >/dev/null 2>&1

