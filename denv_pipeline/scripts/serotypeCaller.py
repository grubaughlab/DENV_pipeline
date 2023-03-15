#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO
from Bio import SeqRecord

#all the names wnat to come out as Yale-XXXX/Yale-XXXX.DENVX.20.cons.fa at this point
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment")
    parser.add_argument("--bed-file", dest="bed_file")

    args = parser.parse_args()

    minCov=50
    amb_list = {"n", "-", "N"}

    #curent ID is formatted Yale-DN00931/Yale-DN00931.DENV1.20.cons.fa 
    if os.path.exists(args.alignment): #when would .cons.fa not exist? is it only if it breaks?

        for sequence in SeqIO.parse(args.alignment, 'fasta'):
            if sequence.seq:
                seqLen=len(sequence.seq)
                seqLenNoAmb = len([i for i in sequence.seq if i not in amb_list])
                percCov=round(seqLenNoAmb/seqLen*100,2)

                if args.bed_file:
                    if os.path.exists(args.bed_file):
                        with open(args.bed_file) as f:
                            for l in f:
                                trimPos = [int(i) for i in l.strip("\n").split("\t")]
                        seqTrim=sequence.seq[int(trimPos[0])-1:int(trimPos[1])]
                        seqLenTrim=len(seqTrim)
                        seqLenNoAmbTrim=len([i for i in seqTrim if i not in amb_list])
                        percCovTrim=round(seqLenNoAmbTrim/seqLenTrim*100,2)
                        
                        if percCovTrim>=minCov:
                            print(f"{sequence.id}\t{seqLen}\t{seqLenNoAmb}\t{percCov}\t{percCovTrim}")
                        else:
                            print(f'{sequence.id}\tNA\tNA\t{seqLen}\t{seqLenNoAmb}\t{percCov}\t{percCovTrim}') 

                    
                    new_header = sequence.id.split("/")[0]
                    seqTrim = SeqRecord.SeqRecord(seqTrim)
                    seqTrim.id = new_header
                    seqTrim.name = new_header
                    seqTrim.description = new_header
                    with open(args.alignment.replace(".out.aln",".out.trim.aln"), 'w') as new_file:
                        SeqIO.write(seqTrim, new_file, 'fasta')

                else:
                    if percCov>=minCov:
                        print(f"{sequence.id}\t{seqLen}\t{seqLenNoAmb}\t{percCov}\tNA")
                    else:
                        print(f'{sequence.id}\tNA\tNA\t{seqLen}\t{seqLenNoAmb}\t{percCov}\tNA')  

            else:
                print(f"{sequence.id}/Unknown\tNA\tNA\tNA\t0\t0")
    else:
        filename = f"{'.'.join(args.alignment.split('.')[0:2])}.cons.fa"
        print(f"{filename}/Unknown\tNA\tNA\tNA\t0\t0")

if __name__=="__main__":
    main()


