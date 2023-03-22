#!/usr/bin/env bash

outdir=$1
#bash
echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > DENV.serotype.calls.tsv
cat ${outdir}/tmp.*.serotype.calls.*.txt >> ${outdir}/DENV.serotype.calls.tsv

#without the depth, these are the same thing, so just renaming this - bash
#not sure star to star is gonna work, might need to do this better
mv ${outdir}/*.*.serotype.txt > ${outdir}/*.serotype.calls.txt

#move to python
echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > ${outdir}/DENV.serotype.calls.mincov50.final.tsv
cat ${outdir}/DENV.serotype.calls.tsv | awk '{ if( $8>=50 ){ print } }' >>  ${outdir}/DENV.serotype.calls.mincov50.final.tsv

#bash
echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" >  ${outdir}/summary.all.samples.tsv
cat ${outdir}/*.serotype.calls.txt >> ${outdir}/summary.all.samples.tsv

#from here on, python
echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > ${outdir}/DENV.top.serotype.calls.all.samples.tsv;

#this sorts the serotype file per sample by coverage and takes the top one, writes it to file
ls ${outdir}/*.serotype.calls.txt | while read i; 
    do 
        cat $i | sort -k8 -n -r | head -1 >> ${outdir}/DENV.top.serotype.calls.all.samples.tsv; 
    done


#copies
cp DENV.top.serotype.calls.all.samples.tsv FINAL;

