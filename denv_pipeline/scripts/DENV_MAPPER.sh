#!/usr/bin/env bash
fname=$1
read1=$2
read2=$3
primer_dir=$4
serotype_caller=$5
depth=$6
outdir=$7
log=$8

#sort out output files - probably mostly change names so they make sense
while IFS= read -r virustype || [[ -n "$virustype" ]]; do 

    fasta=${primer_dir}${virustype}.fasta
    bed=${primer_dir}${virustype}.bed
    trimbed=${primer_dir}${virustype}.trim.bed

    if ! [ -f ${primer_dir}${virustype}.fasta.ann ]; then
        echo "----->>>> making index files for ${virustype}"
        bwa index ${fasta}
    fi

    echo "----->>>>>Mapping reads against serotype "${virustype}" reference sequence"
    bwa mem -v 1 -t 16 ${fasta} $read1 $read2 | samtools view -bS -F 4 -F 2048 | samtools sort -o ${outdir}/${fname%.*}.${virustype}.bam >> ${log} 2>&1

    echo "----->>>>>Trimming bam file"
    ivar trim -e -i ${outdir}/${fname%.*}.${virustype}.bam -b ${bed} -p ${outdir}/${fname%.*}.${virustype}.trimmed.bam >> ${log} 2>&1

    if ! [ -s  ${outdir}/${fname%.*}.${virustype}.trimmed.bam ]; then
        echo "no trimmed bam file found, likely because no reads mapped successfully, exiting shell script"
        continue
    fi

    echo "----->>>>>Sorting bam file"
    samtools sort ${outdir}/${fname%.*}.${virustype}.trimmed.bam -o ${outdir}/${fname%.*}.${virustype}.sort.bam >> ${log} 2>&1

    echo "----->>>>>Indexing bam file"
    samtools index ${outdir}/${fname%.*}.${virustype}.sort.bam >> ${log} 2>&1

    echo "----->>>>>Generating consensus sequence"
    samtools mpileup -aa --reference ${fasta} -A -d 10000 -Q 0 ${outdir}/${fname%.*}.${virustype}.sort.bam | ivar consensus -t 0.75 -m ${depth} -p ${outdir}/${fname%.*}.${virustype}.${depth}.cons -i ${fname%.*}"/"${fname%.*}.${virustype}.${depth}.cons.fa >> ${log} 2>&1
    
    echo "----->>>>>Aligning consensus cps sequence against the reference serotype "${virustype}" cps sequence"
    nextalign run  --reference ${fasta} --output-fasta ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln ${outdir}/${fname%.*}.${virustype}.${depth}.cons.fa >> ${log} 2>&1
    
    if [ -s ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln ]; then
        echo "----->>>>>Aligning with nextalign successful"
    else
        echo "----->>>>>Aligning with mafft (nextalign not successful)" 
        mafft --quiet --6merpair --keeplength  --addfragments ${outdir}/${fname%.*}.${virustype}.${depth}.cons.fa ${fasta} > ${outdir}/ZZ.tmp000.${fname%.*}.${virustype}.${depth}
        grep -A 30000000 `grep ">" ${outdir}/${fname%.*}.${virustype}.${depth}.cons.fa` ${outdir}/ZZ.tmp000.${fname%.*}.${virustype}.${depth} > ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln
    fi

    echo "----->>>>>Calculating percentage coverage against serotype "${virustype}" cps reference sequence"
    if [ -s ${trimbed} ]; then
        python ${serotype_caller} --alignment ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln --bed-file ${trimbed} >> ${outdir}/${fname%.*}_all_virustype_info.txt
    else
        python ${serotype_caller}  --alignment ${outdir}/${fname%.*}.${virustype}.${depth}.out.aln  >> ${outdir}/${fname%.*}_all_virustype_info.txt
    fi

    echo "----->>>>>Identifying variants"
    samtools mpileup -aa --reference ${fasta} -A -d 0 -Q 0 ${outdir}/${fname%.*}.${virustype}.sort.bam | ivar variants -p ${outdir}/${fname%.*}.${virustype}.${depth}.variants -q 20 -t 0.03 -r ${fasta} >> ${log} 2>&1
    
    bedtools genomecov -d -ibam ${outdir}/${fname%.*}.${virustype}.sort.bam > ${outdir}/${fname%.*}.${virustype}.depth.txt; 

done < "${primer_dir}refs.txt"




