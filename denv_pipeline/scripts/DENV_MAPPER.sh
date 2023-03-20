#!/usr/bin/env bash
fname=$1
read1=$2
read2=$3
primer_dir=$4
serotype_caller=$5
outdir=$6

depth=20 #check with chantal if this makes sense for other people

#do proper log files - name them by the sample name put them in a folder.
#currently the log gets written to dsq-jobs-XXXX.out on HPC - either change dsq file or rename jobs file afterwards
#everything needs to be written to outdir
cat ${primer_dir}/DENV.refs.txt | while read denvtype; do 

    fasta=${primer_dir}/${denvtype}.fasta
    bed=${primer_dir}/${denvtype}.bed
    trimbed=${primer_dir}/${denvtype}.trim.bed

    echo "----->>>>>Mapping reads against serotype "${denvtype}" reference sequence"
    bwa mem -v 1 -t 16 ${fasta} $read1 $read2 | samtools view -bS -F 4 -F 2048 | samtools sort -o ${outdir}/${fname%.*}.${denvtype}.bam #> /dev/null 2>&1

    echo "----->>>>>Trimming bam file"
    ivar trim -e -i ${outdir}/${fname%.*}.${denvtype}.bam -b ${bed} -p ${outdir}/${fname%.*}.${denvtype}.trimmed.bam #> /dev/null 2>&1

    echo "----->>>>>Sorting bam file"
    samtools sort ${outdir}/${fname%.*}.${denvtype}.trimmed.bam -o ${outdir}/${fname%.*}.${denvtype}.sort.bam #> /dev/null 2>&1

    echo "----->>>>>Indexing bam file"
    samtools index ${outdir}/${fname%.*}.${denvtype}.sort.bam #> /dev/null 2>&1

#where the loop for depth starts
    echo "----->>>>>Generating consensus sequence"
    samtools mpileup -aa --reference ${fasta} -A -d 10000 -Q 0 ${outdir}/${fname%.*}.${denvtype}.sort.bam | ivar consensus -t 0.75 -m ${depth} -p ${outdir}/${fname%.*}.${denvtype}.${depth}.cons -i ${fname%.*}"/"${fname%.*}.${denvtype}.${depth}.cons".fa" > /dev/null 2>&1
    
    echo "----->>>>>Aligning consensus cps sequence against the reference serotype "${denvtype}" cps sequence"
    nextalign run  --reference ${fasta} --output-fasta ${outdir}/${fname%.*}.${denvtype}.${depth}.out.aln ${outdir}/${fname%.*}.${denvtype}.${depth}.cons.fa
    
    if [ -f ${outdir}/${fname%.*}.${denvtype}.${depth}.out.aln ]; then
        echo "----->>>>>Aligning with nextalign successful"
    else
        echo "----->>>>>Aligning with mafft (nextalign not successful)"
        mafft --quiet --6merpair --maxambiguous 0.99 --addfragments --keeplength --addfragments ${outdir}/${fname%.*}.${denvtype}.${depth}.cons.fa ${ref} > ${outdir}/ZZ.tmp000.${fname%.*}.${denvtype}.${depth}
        grep -A 30000000 `grep ">" ${fname%.*}.${denvtype}.${depth}.cons.fa` ${outdir}/ZZ.tmp000.${fname%.*}.${denvtype}.${depth} > ${outdir}/${fname%.*}.${denvtype}.${depth}.out.aln
    fi

    echo "----->>>>>Calculating percentage coverage against serotype "${denvtype}" cps reference sequence"
    if [ -s ${trimbed} ]; then
        echo "python ${serotype_caller} --alignment ${outdir}/${fname%.*}.${denvtype}.${depth}.out.aln --bed-file ${trimbed} >> ${outdir}/${fname%.*}.${depth}.serotype.txt"
        python ${serotype_caller} --alignment ${outdir}/${fname%.*}.${denvtype}.${depth}.out.aln --bed-file ${trimbed} >> ${outdir}/${fname%.*}.${depth}.serotype.txt
    else
        python ${serotype_caller}  --alignment ${outdir}/${fname%.*}.${dentype}.${depth}.out.aln  >> ${outdir}/${fname%.*}.${depth}.serotype.txt
    fi

    echo "----->>>>>Identifying variants"
    samtools mpileup -aa --reference ${fasta} -A -d 0 -Q 0 ${outdir}/${fname%.*}.${denvtype}.sort.bam | ivar variants -p ${outdir}/${fname%.*}.${denvtype}.${depth}_variants -q 20 -t 0.03 -r ${fasta} 
    
    #pulls variants which are above 20% and below 80% and have "PASS" as TRUE
    awk -F '\t' '(0.2 < $11 && $11 < 0.80) && ($14=="TRUE")' ${outdir}/${fname%.*}.${denvtype}.${depth}_variants.tsv > ${outdir}/${fname%.*}.${denvtype}.${depth}_variants_frequency.tsv
    echo -e ${fname%.*}"\t"`wc -l ${outdir}/${fname%.*}.${denvtype}.${depth}_variants_frequency.tsv | awk '{print $1}'` > ${outdir}/${fname%.*}.${denvtype}.${depth}_variants_frequency_count.txt

    #depth loop finishes here
    
    bedtools genomecov -d -ibam ${outdir}/${fname%.*}.${denvtype}.sort.bam > ${outdir}/${fname%.*}.${denvtype}.depth.txt; 

done


#other depth loop here
cat ${outdir}/${fname%.*}.${depth}.serotype.txt | sort -k8 -n -r | awk '{ if( $8>=50 ){ print } }' >> ${outdir}/tmp.${fname%.*}.serotype.calls.${depth}.txt

cat ${outdir}/tmp.${fname%.*}.serotype.calls.*.txt >> ${outdir}/DENV.serotype.calls.tsv
cat ${outdir}/${fname%.*}.*.serotype.txt > ${outdir}/${fname%.*}.serotype.calls.txt


