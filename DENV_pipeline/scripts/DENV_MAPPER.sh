#!/usr/bin/env bash

read1=$1
read2=$2
cps_refs=$3
primer_bed=$4

fname=`basename ${read1} | tr "._" "\t\t" | awk '{print $1}'`

rm -rf ${fname%.*}.serotype.txt > /dev/null 2>&1
rm -rf tmp.${fname%.*}.serotype.calls.*.txt

cat ${cps_refs} | while read ref; do 
    cps=`basename ${ref} | tr "\_\t." "\t\t" | awk '{print $1}'`
    cps1=`dirname ${ref}`"/"`basename ${ref} | tr "\_\t\." "\t\t\t" | awk '{print $1}'`

    echo "----->>>>>Mapping reads against serotype "${cps}" reference sequence"
    bwa mem -v 1 -t 16 ${ref} $read1 $read2 | samtools view -bS -F 4 -F 2048 | samtools sort -o ${fname%.*}.${cps}.bam > /dev/null 2>&1

    echo "----->>>>>Trimming bam file"
    ivar trim -e -i ${fname%.*}.${cps}.bam -b ${cps1}.bed -p ${fname%.*}.${cps}.trimmed.bam > /dev/null 2>&1

    echo "----->>>>>Sorting bam file"
    samtools sort ${fname%.*}.${cps}.trimmed.bam -o ${fname%.*}.${cps}.sort.bam > /dev/null 2>&1

    echo "----->>>>>Indexing bam file"
    samtools index ${fname%.*}.${cps}.sort.bam > /dev/null 2>&1

    for depth in 20
    do
        echo "----->>>>>Generating consensus sequence"
        samtools mpileup -aa --reference ${ref} -A -d 10000 -Q 0 ${fname%.*}.${cps}.sort.bam | ivar consensus -t 0.75 -m ${depth} -p ${fname%.*}.${cps}.${depth}.cons -i ${fname%.*}"/"${fname%.*}.${cps}.${depth}.cons".fa/"${cps}"/"${ref} > /dev/null 2>&1

        rm -rf ${fname%.*}.${cps}.out.aln  > /dev/null 2>&1
        echo "----->>>>>Aligning consensus cps sequence against the reference serotype "${cps}" cps sequence"

        echo "nextalign --output-fasta ${fname%.*}.${cps}.${depth}.out.aln --reference ${ref} --sequences ${fname%.*}.${cps}.${depth}.cons.fa"

        nextalign --output-fasta ${fname%.*}.${cps}.${depth}.out.aln --reference ${ref} --sequences ${fname%.*}.${cps}.${depth}.cons.fa > /dev/null 2>&1

        if [ -f ${fname%.*}.${cps}.out.aln ]; then
            echo "----->>>>>Aligning with nextalign successful"
        else
            echo "----->>>>>Aligning with mafft (nextalign not successful)"
            mafft --quiet --6merpair --maxambiguous 0.99 --addfragments --keeplength --addfragments ${fname%.*}.${cps}.${depth}.cons.fa ${ref} > ZZ.tmp000.${fname%.*}.${cps}.${depth}
            grep -A 30000000 `grep ">" ${fname%.*}.${cps}.${depth}.cons.fa` ZZ.tmp000.${fname%.*}.${cps}.${depth} > ${fname%.*}.${cps}.${depth}.out.aln
        fi
        
        rm -rf ZZ.tmp000.${fname%.*}.${cps}.${depth}

        echo "----->>>>>Calculating percentage coverage against serotype "${cps}" cps reference sequence"


        if [ -s ${cps}.trim.bed ]; then
            python /gpfs/ycga/project/grubaugh/shared/DENVSEQ/SCRIPT/serotypeCaller.py ${fname%.*}.${cps}.${depth}.out.aln ${cps1}.trim.bed | tr "/" "\t" | awk -v depth="$depth" '{print $1"\t"$2"\t"depth"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' >> ${fname%.*}.${depth}.serotype.txt
        else
            python /gpfs/ycga/project/grubaugh/shared/DENVSEQ/SCRIPT/serotypeCaller.py ${fname%.*}.${cps}.${depth}.out.aln | tr "/" "\t" | awk -v depth="$depth" '{print $1"\t"$2"\t"depth"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' >> ${fname%.*}.${depth}.serotype.txt
        fi

        echo "----->>>>>Identifying variants"
        samtools mpileup -aa --reference ${ref} -A -d 0 -Q 0 ${fname%.*}.${cps}.sort.bam | ivar variants -p ${fname%.*}.${cps}.${depth}_variants -q 20 -t 0.03 -r ${ref} 

        awk -F '\t' '(0.2 < $11 && $11 < 0.80) && ($14=="TRUE")' ${fname%.*}.${cps}.${depth}_variants.tsv > ${fname%.*}.${cps}.${depth}_variants_frequency.tsv
        echo -e ${fname%.*}"\t"`wc -l ${fname%.*}.${cps}.${depth}_variants_frequency.tsv | awk '{print $1}'` > ${fname%.*}.${cps}.${depth}_variants_frequency_count.txt
    
        rm -rf ${fname%.*}.${cps}.cons.qual.txt ${fname%.*}.${cps}.bam ${fname_path} ${fname%.*}.${cps}.sort.bam*bai ${fname%.*}.${cps}.trimmed.bam
    done
    
    bedtools genomecov -d -ibam ${fname%.*}.${cps}.sort.bam > ${fname%.*}.${cps}.depth.txt; 
done

rm -rf tmp.${fname%.*}.serotype.calls.*.txt > /dev/null 2>&1

for depth in 20
do
    cat ${fname%.*}.${depth}.serotype.txt | sort -k8 -n -r | awk '{ if( $8>=50 ){ print } }' >> tmp.${fname%.*}.serotype.calls.${depth}.txt
done

cat tmp.${fname%.*}.serotype.calls.*.txt >> DENV.serotype.calls.tsv
rm -rf tmp.${fname%.*}.serotype.calls.*.txt > /dev/null 2>&1

cat ${fname%.*}.*.serotype.txt > ${fname%.*}.serotype.calls.txt
rm -rf ${fname%.*}.*.serotype.txt

