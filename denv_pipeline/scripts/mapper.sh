#!/usr/bin/env bash
fname=$1
read1=$2
read2=$3
primer_dir=$4
serotype_caller=$5
depth=$6
tempdir=$7
log=$8

while IFS= read -r virustype || [[ -n "$virustype" ]]; do 

    fasta=${primer_dir}/${virustype}.fasta
    bed=${primer_dir}/${virustype}.bed
    trimbed=${primer_dir}/${virustype}.trim.bed
    consensus_name=${fname}.${virustype}

    if ! [ -f ${primer_dir}/${virustype}.fasta.ann ]; then
        echo "----->>>> making index files for ${virustype}"
        bwa index ${fasta} >> ${log} 2>&1
    fi

    echo "----->>>>>Mapping reads against serotype "${virustype}" reference sequence"
    bwa mem -v 1 -t 16 ${fasta} $read1 $read2 > bwa_mem_output.txt >> ${log} 2>&1
    samtools view -bS -F 4 -F 2048 bwa_mem_output.txt | samtools sort -o ${tempdir}/${fname}.${virustype}.bam >> ${log} 2>&1

    echo "----->>>>>Trimming bam file"
    ivar trim -e -i ${tempdir}/${fname}.${virustype}.bam -b ${bed} -p ${tempdir}/${fname}.${virustype}.trimmed.bam >> ${log} 2>&1

    if ! [ -s  ${tempdir}/${fname%.*}.${virustype}.trimmed.bam ]; then
        echo "no trimmed bam file found, likely because no reads mapped successfully, exiting shell script"
        python -c "import denv_pipeline.scripts.make_summary_files; denv_pipeline.scripts.make_summary_files.make_empty_file('${tempdir}', '${depth}', '${fname}', '${virustype}')"
        continue
    fi

    echo "----->>>>>Sorting bam file"
    samtools sort ${tempdir}/${fname}.${virustype}.trimmed.bam -o ${tempdir}/${fname}.${virustype}.sort.bam >> ${log} 2>&1

    echo "----->>>>>Indexing bam file"
    samtools index ${tempdir}/${fname}.${virustype}.sort.bam >> ${log} 2>&1

    echo "----->>>>>Generating consensus sequence"
    samtools mpileup -aa --reference ${fasta} -A -d 10000 -Q 0 ${tempdir}/${fname}.${virustype}.sort.bam | ivar consensus -t 0.75 -m ${depth} -p ${tempdir}/${fname}.${virustype}.${depth}.cons -i ${consensus_name} >> ${log} 2>&1
    
    echo "----->>>>>Aligning consensus cps sequence against the reference serotype "${virustype}" cps sequence"
    nextalign run  --reference ${fasta} --output-fasta ${tempdir}/${fname}.${virustype}.${depth}.out.aln ${tempdir}/${fname}.${virustype}.${depth}.cons.fa >> ${log} 2>&1
    
    if [ -s ${tempdir}/${fname%.*}.${virustype}.${depth}.out.aln ]; then
        echo "----->>>>>Aligning with nextalign successful"
    else
        echo "----->>>>>Aligning with mafft (nextalign not successful)" 
        mafft --quiet --6merpair --keeplength  --addfragments ${tempdir}/${fname}.${virustype}.${depth}.cons.fa ${fasta} > ${tempdir}/${fname}.${virustype}.${depth}.alignment_intermediate.fasta
        grep -A 30000000 ">${consensus_name}" ${tempdir}/${fname}.${virustype}.${depth}.alignment_intermediate.fasta > ${tempdir}/${fname}.${virustype}.${depth}.out.aln
    fi

    echo "----->>>>>Calculating percentage coverage against serotype "${virustype}" cps reference sequence"
    if [ -s ${trimbed} ]; then
        python ${serotype_caller} --alignment ${tempdir}/${fname}.${virustype}.${depth}.out.aln --bed-file ${trimbed} --outfile ${tempdir}/${fname}_all_virustype_info.txt
    else
        python ${serotype_caller}  --alignment ${tempdir}/${fname}.${virustype}.${depth}.out.aln --outfile ${tempdir}/${fname}_all_virustype_info.txt
    fi

    echo "----->>>>>Identifying variants"
    samtools mpileup -aa --reference ${fasta} -A -d 0 -Q 0 ${tempdir}/${fname}.${virustype}.sort.bam | ivar variants -p ${tempdir}/${fname}.${virustype}.${depth}.variants -q 20 -t 0.03 -r ${fasta} >> ${log} 2>&1
    
    echo "----->>>>>>Getting depths"
    bedtools genomecov -d -ibam ${tempdir}/${fname}.${virustype}.sort.bam > ${tempdir}/${fname}.${virustype}.depth.txt; 

    echo "--->>>>> Finished"

done < "${primer_dir}/refs.txt"
