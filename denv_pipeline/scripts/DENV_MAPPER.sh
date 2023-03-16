#!/usr/bin/env bash
fname=$1
read1=$2
read2=$3
primer_dir=$4
primer_bed=$5


cat ${primer_dir}DENV.refs.txt | while read denvtype; do 

    depth=20
    fasta=${primer_dir}${denvtype}.fasta
    bed=${primer_dir}${denvtype}.bed

    echo "----->>>>>Mapping reads against serotype "${denvtype}" reference sequence"
    bwa mem -v 1 -t 16 ${fasta} $read1 $read2 | samtools view -bS -F 4 -F 2048 | samtools sort -o ${fname%.*}.${denvtype}.bam #> /dev/null 2>&1

    echo "----->>>>>Trimming bam file"
    ivar trim -e -i ${fname%.*}.${denvtype}.bam -b ${bed} -p ${fname%.*}.${denvtype}.trimmed.bam #> /dev/null 2>&1

    echo "----->>>>>Sorting bam file"
    samtools sort ${fname%.*}.${denvtype}.trimmed.bam -o ${fname%.*}.${denvtype}.sort.bam #> /dev/null 2>&1

    echo "----->>>>>Indexing bam file"
    samtools index ${fname%.*}.${denvtype}.sort.bam #> /dev/null 2>&1

#where the loop for depth starts
    echo "----->>>>>Generating consensus sequence"
    samtools mpileup -aa --reference ${fasta} -A -d 10000 -Q 0 ${fname%.*}.${denvtype}.sort.bam | ivar consensus -t 0.75 -m ${depth} -p ${fname%.*}.${denvtype}.${depth}.cons -i ${fname%.*}"/"${fname%.*}.${denvtype}.${depth}.cons".fa" > /dev/null 2>&1
    
    echo "----->>>>>Aligning consensus cps sequence against the reference serotype "${denvtype}" cps sequence"
    nextalign run  --reference ${fasta} --output-fasta ${fname%.*}.${denvtype}.${depth}.out.aln ${fname%.*}.${denvtype}.${depth}.cons.fa
    
    if [ -f ${fname%.*}.${denvtype}.${depth}.out.aln ]; then
        echo "----->>>>>Aligning with nextalign successful"
    else
        echo "----->>>>>Aligning with mafft (nextalign not successful)"
        mafft --quiet --6merpair --maxambiguous 0.99 --addfragments --keeplength --addfragments ${fname%.*}.${denvtype}.${depth}.cons.fa ${ref} > ZZ.tmp000.${fname%.*}.${denvtype}.${depth}
        grep -A 30000000 `grep ">" ${fname%.*}.${denvtype}.${depth}.cons.fa` ZZ.tmp000.${fname%.*}.${denvtype}.${depth} > ${fname%.*}.${denvtype}.${depth}.out.aln
    fi

    rm ZZ.tmp000.${fname%.*}.${denvtype}.${depth} 2>&1

    echo "----->>>>>Calculating percentage coverage against serotype "${denvtype}" cps reference sequence"
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
#depth loop finishes here
    
    bedtools genomecov -d -ibam ${fname%.*}.${cps}.sort.bam > ${fname%.*}.${cps}.depth.txt; 


rm -rf tmp.${fname%.*}.serotype.calls.*.txt > /dev/null 2>&1

for depth in 20
do
    cat ${fname%.*}.${depth}.serotype.txt | sort -k8 -n -r | awk '{ if( $8>=50 ){ print } }' >> tmp.${fname%.*}.serotype.calls.${depth}.txt
done

cat tmp.${fname%.*}.serotype.calls.*.txt >> DENV.serotype.calls.tsv
rm -rf tmp.${fname%.*}.serotype.calls.*.txt > /dev/null 2>&1

cat ${fname%.*}.*.serotype.txt > ${fname%.*}.serotype.calls.txt
rm -rf ${fname%.*}.*.serotype.txt

