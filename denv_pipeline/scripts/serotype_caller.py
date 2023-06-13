#!/usr/bin/env python

import os
import sys
import argparse
import csv
from Bio import SeqIO
from Bio import SeqRecord

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment")
    parser.add_argument("--bed-file", dest="bed_file")
    parser.add_argument("--outfile")

    args = parser.parse_args()
    headers = ["sample_id","consensus_sequence_file","depth","serotype","reference_sequence_name","reference_sequence_length","number_aligned_bases","coverage_untrimmed","coverage_trimmed"]

    if not os.path.exists(args.outfile):
        with open(args.outfile, 'w') as fw:
            writer = csv.DictWriter(fw, delimiter = "\t", fieldnames=headers)
            writer.writeheader()
            write_dict = populate_line(args)

            writer.writerow(write_dict)
        
    else:        
        with open(args.outfile, 'a') as fw:
            writer = csv.DictWriter(fw, delimiter = "\t", fieldnames=headers)
            write_dict = populate_line(args)
            writer.writerow(write_dict)


def populate_line(args):

    min_coverage=50
    amb_list = {"n", "-", "N"}

    name_elements = args.alignment.split("/")[-1].split(".")
    serotype = name_elements[1]
    ref_sequence = serotype
    
    write_dict = {}
    write_dict["sample_id"] = name_elements[0]
    write_dict["consensus_sequence_file"] = f'{".".join(name_elements[0:3])}.cons.fa'
    write_dict["depth"] = name_elements[2]
    write_dict["reference_sequence_name"] = ref_sequence

    if os.path.exists(args.alignment) and os.path.getsize(args.alignment) > 0: 
        for sequence in SeqIO.parse(args.alignment, 'fasta'):
            if sequence.seq:
                seq_len=len(sequence.seq)
                seq_len_no_amb = len([i for i in sequence.seq if i not in amb_list])
                perc_cov=round(seq_len_no_amb/seq_len*100,2)

                write_dict["reference_sequence_length"] = seq_len
                write_dict["number_aligned_bases"] = seq_len_no_amb
                write_dict["coverage_untrimmed"] = perc_cov

                if args.bed_file:
                    if os.path.exists(args.bed_file):
                        with open(args.bed_file) as f:
                            for l in f:
                                trim_pos = [int(i) for i in l.strip("\n").split("\t")]
                        seq_trim = sequence.seq[int(trim_pos[0])-1:int(trim_pos[1])]
                        seq_len_trim = len(seq_trim)
                        seq_len_no_amb_trim=len([i for i in seq_trim if i not in amb_list])
                        perc_cov_trim=round(seq_len_no_amb_trim/seq_len_trim*100,2)

                        write_dict["coverage_trimmed"] = perc_cov_trim

                        if perc_cov_trim>=min_coverage:
                            write_dict["serotype"] = serotype
                        else:
                            write_dict["serotype"] = "NA"
                    
                        seq_trim = SeqRecord.SeqRecord(seq_trim)
                        with open(args.alignment.replace(".out.aln",".out.trim.aln"), 'w') as new_file: 
                            SeqIO.write(seq_trim, new_file, 'fasta')
                    else:
                        sys.stderr.write(f"Bed file {args.bed_file} not found")
                        with open(args.alignment.replace(".out.aln",".out.trim.aln"), 'w') as new_file:
                            pass


                else:
                    write_dict["coverage_trimmed"] = 0
                    if perc_cov>=min_coverage:
                        write_dict["serotype"] = serotype
                    else:
                        write_dict["serotype"] = "NA"
                    with open(args.alignment.replace(".out.aln",".out.trim.aln"), 'w') as new_file:
                        pass

            else:
                write_dict["serotype"] = "NA"
                write_dict["reference_sequence_length"] = "NA"
                write_dict["number_aligned_bases"] = 0
                write_dict["coverage_untrimmed"] = 0
                write_dict["coverage_trimmed"] = 0
                
    else:
        write_dict["serotype"] = "NA"
        write_dict["reference_sequence_length"] = "NA"
        write_dict["number_aligned_bases"] = 0
        write_dict["coverage_untrimmed"] = 0
        write_dict["coverage_trimmed"] = 0


    return write_dict

if __name__=="__main__":
    main()
