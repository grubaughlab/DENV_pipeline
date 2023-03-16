#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO
from Bio import SeqRecord

#all the names wnat to come out as Yale-XXXX/Yale-XXXX.DENVX.20.cons.fa at this point
#make everything snakecase?
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment")
    parser.add_argument("--bed-file", dest="bed_file")

    args = parser.parse_args()

    min_coverage=50
    amb_list = {"n", "-", "N"}

    #curent ID is formatted Yale-DN00931/Yale-DN00931.DENV1.20.cons.fa 
    if os.path.exists(args.alignment): #when would .cons.fa not exist? is it only if it breaks?

        for sequence in SeqIO.parse(args.alignment, 'fasta'):
            if sequence.seq:
                seq_len=len(sequence.seq)
                seq_len_no_amb = len([i for i in sequence.seq if i not in amb_list])
                perc_cov=round(seq_len_no_amb/seq_len*100,2)

                if args.bed_file:
                    if os.path.exists(args.bed_file):
                        with open(args.bed_file) as f:
                            for l in f:
                                trim_pos = [int(i) for i in l.strip("\n").split("\t")]
                        seq_trim=sequence.seq[int(trim_pos[0])-1:int(trim_pos[1])]
                        seq_len_trim=len(seq_trim)
                        seq_len_no_amb_trim=len([i for i in seq_trim if i not in amb_list])
                        perc_cov_trim=round(seq_len_no_amb_trim/seq_len_trim*100,2)
                        
                        if perc_cov_trim>=min_coverage:
                            print(f"{sequence.id}\t{seq_len}\t{seq_len_no_amb}\t{perc_cov}\t{perc_cov_trim}")
                        else:
                            print(f'{sequence.id}\tNA\tNA\t{seq_len}\t{seq_len_no_amb}\t{perc_cov}\t{perc_cov_trim}') 

                    
                    new_header = sequence.id.split("/")[0]
                    seq_trim = SeqRecord.SeqRecord(seq_trim)
                    seq_trim.id = new_header
                    seq_trim.name = new_header
                    seq_trim.description = new_header
                    with open(args.alignment.replace(".out.aln",".out.trim.aln"), 'w') as new_file:
                        SeqIO.write(seq_trim, new_file, 'fasta')

                else:
                    if perc_cov>=min_coverage:
                        print(f"{sequence.id}\t{seq_len}\t{seq_len_no_amb}\t{perc_cov}\tNA")
                    else:
                        print(f'{sequence.id}\tNA\tNA\t{seq_len}\t{seq_len_no_amb}\t{perc_cov}\tNA')  

            else:
                print(f"{sequence.id}/Unknown\tNA\tNA\tNA\t0\t0")
    else:
        filename = f"{'.'.join(args.alignment.split('.')[0:2])}.cons.fa"
        print(f"{filename}/Unknown\tNA\tNA\tNA\t0\t0")

if __name__=="__main__":
    main()


