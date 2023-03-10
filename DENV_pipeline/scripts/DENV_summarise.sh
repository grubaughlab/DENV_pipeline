#!/usr/bin/env bash

rm -rf FINAL > /dev/null 2>&1

echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > DENV.serotype.calls.final.tsv
cat DENV.serotype.calls.tsv >> DENV.serotype.calls.final.tsv

echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > DENV.serotype.calls.mincov50.final.tsv
cat DENV.serotype.calls.tsv | awk '{ if( $8>=50 ){ print } }' >> DENV.serotype.calls.mincov50.final.tsv

echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > summary.all.samples.tsv
cat *.serotype.calls.txt >> summary.all.samples.tsv

echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > DENV.top.serotype.calls.all.samples.tsv;

ls *.serotype.calls.txt | while read i; 
    do 
        cat $i | sort -k8 -n -r | head -1 >> DENV.top.serotype.calls.all.samples.tsv; 
    done

mkdir -p FINAL/BAM FINAL/VARIANTS FINAL/DEPTH FINAL/ALIGNMENT FINAL/CONSENSUS;

cp DENV.serotype.calls.final.tsv FINAL;
cp summary.all.samples.tsv FINAL;

mkdir -p FINAL/VARIANTS FINAL/DEPTH FINAL/ALIGNMENT FINAL/CONSENSUS;

grep -v ConsSequence DENV.serotype.calls.mincov50.final.tsv | while read i; 
    do 
        cp `echo -e $i | awk '{print $1"."$4".sort.bam"}'` FINAL/BAM/; 
    done

grep -v ConsSequence DENV.serotype.calls.mincov50.final.tsv | cut -f2 | while read i; 
    do 
        cat ${i} | tr "/" "\t" | awk '{print $1}' > FINAL/CONSENSUS/${i}; 
    done

grep -v ConsSequence DENV.serotype.calls.mincov50.final.tsv | cut -f2 | while read i; 
    do 
        NN=`echo $i | cut -f2 | tr "." "\t" | awk '{print $2}'`; 
        j=${i%.*}; 
        r=${j%.*}; 
        cat ${r}.out.trim.aln >> FINAL/ALIGNMENT/${NN}.trim.aln; 
    done

grep -v ConsSequence DENV.serotype.calls.mincov50.final.tsv | cut -f2 | while read i; 
    do 
        NN=`echo $i | cut -f2 | tr "." "\t" | awk '{print $2}'`; 
        j=${i%.*}; r=${j%.*}; 
        cat ${r}.out.aln >> FINAL/ALIGNMENT/${NN}.untrim.aln; 
    done

grep -v ConsSequence DENV.serotype.calls.mincov50.final.tsv | cut -f2 | while read i; 
    do  
        j=${i%.*}; 
        r=${j%.*};
        cp ${r%.*}.depth.txt FINAL/DEPTH/; 
    done

grep -v ConsSequence DENV.serotype.calls.mincov50.final.tsv | cut -f2 | while read i; 
    do  
        j=${i%.*}; 
        r=${j%.*}; 
        cp ${r}_variants_frequency.tsv FINAL/VARIANTS; 
    done

grep -v ConsSequence DENV.serotype.calls.mincov50.final.tsv | cut -f2 | while read i; 
    do  
        j=${i%.*}; 
        r=${j%.*}; 
        cp ${r}_variants_frequency_count.txt FINAL/VARIANTS; 
    done

cp DENV.top.serotype.calls.all.samples.tsv FINAL;

cat *variants_frequency_count.txt > Variants_summary.tsv;
cp Variants_summary.tsv FINAL/VARIANTS; 
