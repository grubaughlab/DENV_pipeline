#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO
from Bio import SeqRecord

#all the names want to come out as Yale-XXXX/Yale-XXXX.DENVX.20.cons.fa at this point
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment")
    parser.add_argument("--bed-file", dest="bed_file")

    args = parser.parse_args()

    min_coverage=50
    amb_list = {"n", "-", "N"}

    headers = ["SampleID", "ConsSequence", "Depth", "Serotype", "RefSerotypeSequence", "RefSeqLength", "AlignedBases", "CoverageUntrimmed", "CoverageTrimmed"]

    name_elements = args.alignment.split("/")[-1].split(".")
    serotype = name_elements[1]
    ref_sequence = f'{serotype}.fasta'
    
    write_dict = {}
    write_dict["SampleID"] = name_elements[0]
    write_dict["ConsSequence"] = f'{".".join(name_elements[0:3])}.cons.fa'
    write_dict["Depth"] = name_elements[2]

    #curent ID is formatted Yale-DN00931/Yale-DN00931.DENV1.20.cons.fa 
    if os.path.exists(args.alignment): #when would .cons.fa not exist? is it only if it breaks? Or no reads?

        for sequence in SeqIO.parse(args.alignment, 'fasta'):
            if sequence.seq:
                seq_len=len(sequence.seq)
                seq_len_no_amb = len([i for i in sequence.seq if i not in amb_list])
                perc_cov=round(seq_len_no_amb/seq_len*100,2)

                write_dict["RefSeqLength"] = seq_len
                write_dict["AlignedBases"] = seq_len_no_amb
                write_dict["CoverageUntrimmed"] = perc_cov

                if args.bed_file:
                    if os.path.exists(args.bed_file):
                        with open(args.bed_file) as f:
                            for l in f:
                                trim_pos = [int(i) for i in l.strip("\n").split("\t")]
                        seq_trim = sequence.seq[int(trim_pos[0])-1:int(trim_pos[1])]
                        seq_len_trim = len(seq_trim)
                        seq_len_no_amb_trim=len([i for i in seq_trim if i not in amb_list])
                        perc_cov_trim=round(seq_len_no_amb_trim/seq_len_trim*100,2)

                        write_dict["CoverageTrimmed"] = perc_cov_trim

                        if perc_cov_trim>=min_coverage:
                            write_dict["Serotype"] = serotype
                            write_dict["RefSerotypeSequence"] = ref_sequence
                        else:
                            write_dict["Serotype"] = "NA"
                            write_dict["RefSerotypeSequence"] = "NA"
                    
                        new_header = sequence.id.split("/")[0]
                        seq_trim = SeqRecord.SeqRecord(seq_trim)
                        seq_trim.id = new_header
                        seq_trim.name = new_header
                        seq_trim.description = new_header
                        with open(args.alignment.replace(".out.aln",".out.trim.aln"), 'w') as new_file:
                            SeqIO.write(seq_trim, new_file, 'fasta')
                    else:
                        sys.stderr.write(f"Bed file {args.bed_file} not found")

                else:
                    write_dict["CoverageTrimmed"] = "NA"
                    if perc_cov>=min_coverage:
                        write_dict["Serotype"] = serotype
                        write_dict["RefSerotypeSequence"] = ref_sequence
                    else:
                        write_dict["Serotype"] = "NA"
                        write_dict["RefSerotypeSequence"] = "NA"

            else:
                write_dict["Serotype"] = "Unknown"
                write_dict["RefSerotypeSequence"] = "NA"
                write_dict["RefSeqLength"] = "NA"
                write_dict["AlignedBases"] = "NA"
                write_dict["CoverageUntrimmed"] = 0
                write_dict["CoverageTrimmed"] = 0
    else:
        write_dict["Serotype"] = "Unknown"
        write_dict["RefSerotypeSequence"] = "NA"
        write_dict["RefSeqLength"] = "NA"
        write_dict["AlignedBases"] = "NA"
        write_dict["CoverageUntrimmed"] = 0
        write_dict["CoverageTrimmed"] = 0

    final = []
    for col in headers:
        final.append(str(write_dict[col]))
    
    print("\t".join(final))

if __name__=="__main__":
    main()
