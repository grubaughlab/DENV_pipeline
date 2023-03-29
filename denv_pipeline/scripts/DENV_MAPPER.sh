#!/usr/bin/env bash
fname=$1
read1=$2
read2=$3
primer_dir=$4
serotype_caller=$5
depth=$6
outdir=$7

#sort out output files - probably mostly change names so they make sense
#do proper log files - name them by the sample name put them in a folder.
#currently the log gets written to dsq-jobs-XXXX.out on HPC - either change dsq file or rename jobs file afterwards
cat ${primer_dir}refs.txt | while read virustype; do 

    fasta=${primer_dir}${virustype}.fasta
    bed=${primer_dir}${virustype}.bed
    trimbed=${primer_dir}${virustype}.trim.bed

    echo "----->>>>>Mapping reads against serotype "${virustype}" reference sequence"
    bwa mem -v 1 -t 16 ${fasta} $read1 $read2 | samtools view -bS -F 4 -F 2048 | samtools sort -o ${outdir}/${fname%.*}.${virustype}.bam > /dev/null 2>&1

    echo "----->>>>>Trimming bam file"
    ivar trim -e -i ${outdir}/${fname%.*}.${virustype}.bam -b ${bed} -p ${outdir}/${fname%.*}.${virustype}.trimmed.bam > /dev/null 2>&1

    echo "----->>>>>Sorting bam file"
    samtools sort ${outdir}/${fname%.*}.${virustype}.trimmed.bam -o ${outdir}/${fname%.*}.${virustype}.sort.bam > /dev/null 2>&1

    echo "----->>>>>Indexing bam file"
    samtools index ${outdir}/${fname%.*}.${virustype}.sort.bam > /dev/null 2>&1

#where the loop for depth starts
    echo "----->>>>>Generating consensus sequence"
    samtools mpileup -aa --reference ${fasta} -A -d 10000 -Q 0 ${outdir}/${fname%.*}.${virustype}.sort.bam | ivar consensus -t 0.75 -m ${depth} -p ${outdir}/${fname%.*}.${virustype}.${depth}.cons -i ${fname%.*}"/"${fname%.*}.${virustype}.${depth}.cons".fa" > /dev/null 2>&1
    
    echo "----->>>>>Aligning consensus cps sequence against the reference serotype "${virustype}" cps sequence"
    nextalign run  --reference ${fasta} --output-fasta ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln ${outdir}/${fname%.*}.${virustype}.${depth}.cons.fa
    
    if [ -f ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln ]; then
        echo "----->>>>>Aligning with nextalign successful"
    else
        echo "----->>>>>Aligning with mafft (nextalign not successful)"
        mafft --quiet --6merpair --maxambiguous 0.99 --addfragments --keeplength --addfragments ${outdir}/${fname%.*}.${virustype}.${depth}.cons.fa ${ref} > ${outdir}/ZZ.tmp000.${fname%.*}.${virustype}.${depth}
        grep -A 30000000 `grep ">" ${fname%.*}.${virustype}.${depth}.cons.fa` ${outdir}/ZZ.tmp000.${fname%.*}.${virustype}.${depth} > ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln
    fi

    #to put depth loop back in here, add depth to output name
    echo "----->>>>>Calculating percentage coverage against serotype "${virustype}" cps reference sequence"
    if [ -s ${trimbed} ]; then
        echo "python ${serotype_caller} --alignment ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln --bed-file ${trimbed} >> ${outdir}/${fname%.*}.serotype.calls.txt"
        python ${serotype_caller} --alignment ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln --bed-file ${trimbed} >> ${outdir}/${fname%.*}.serotype.calls.txt
    else
        python ${serotype_caller}  --alignment ${outdir}/${fname%.*}.${dentype}.${depth}.out.aln  >> ${outdir}/${fname%.*}.serotype.calls.txt
    fi

    echo "----->>>>>Identifying variants"
    samtools mpileup -aa --reference ${fasta} -A -d 0 -Q 0 ${outdir}/${fname%.*}.${virustype}.sort.bam | ivar variants -p ${outdir}/${fname%.*}.${virustype}.${depth}.variants -q 20 -t 0.03 -r ${fasta} 

    #depth loop finishes here
    
    bedtools genomecov -d -ibam ${outdir}/${fname%.*}.${virustype}.sort.bam > ${outdir}/${fname%.*}.${virustype}.depth.txt; 

done


#other depth loop here
#takes the top call if it's over 50% coverage
cat ${outdir}/${fname%.*}.serotype.calls.txt | sort -k8 -n -r | awk '{ if( $8>=50 ){ print } }' >> ${outdir}/tmp.${fname%.*}.serotype.calls.${depth}.txt

#if depth loop goes back in, need this at the end - and change the outputs of the python script to be name.depth.serotype.txt
#cat ${outdir}/${fname%.*}.*.serotype.txt > ${outdir}/${fname%.*}.serotype.calls.txt



