import os
import shutil
import sys
import csv
from collections import defaultdict
import datetime as dt
from unidecode import unidecode
from Bio import SeqIO
import tqdm as tqdm
import argparse

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument("--master-csv", dest="master_csv", help="GLab master csv")
    parser.add_argument("--file-path", dest="file_path", help="Overall path to where genomic data is kept, relative to script")
    parser.add_argument("--genbank-path", dest="genbank_path", help="path to where genbank files will be made, relative to current directory")
    parser.add_argument("--date-submission", dest="date_submission", help="for file name within lab folder")
    parser.add_argument("-h","--help",action="store_true",dest="help")

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)

    master_csv = args.master_csv
    file_path = args.file_path
    genbank_path = args.genbank_path
    date_submission = args.date_submission

    
    #parse the metadata to get everything that's needed
    run_to_seqs = defaultdict(list)
    new_names = {}
    seqs_to_lab = {}
    seqs_to_bam_name = {}

    with open(master_csv) as f:
        data = csv.DictReader(f)
        for l in data:
            if l['Final_Genome'] == "1" and l['GenBank_Accession'] == "":
                
                if "+" in l['NGS_Run_ID']:
                    run = f"{l['NGS_Run_ID'].split('+')[0].rstrip(' ')}CONCAT"
                else:
                    run = l['NGS_Run_ID']
                            
                run_to_seqs[run].append(l['Yale-ID'])
                
                
                if l['Country'] == "USA":
                    location = f'{l["Division_(state)"]}_USA'
                else:
                    location = l["Country"].replace(" ","_")
                
                new_name = f'{l["NGS_Serotype"]}/Human/{location}/{l["Collection_Date"]}/{l["Yale-ID"]}'
                new_names[l['Yale-ID']] = new_name
                
                seqs_to_lab[l['Yale-ID']] = unidecode(l["Lab_Source"].replace(" ","_"))
                
                if "_P" in run:
                    new_id = "-".join(l["Yale-ID"].split("-")[0:2])
                    seqs_to_bam_name[l["Yale-ID"]] = f'{new_id}P.{l["NGS_Serotype"]}.sort.bam'
                else:
                    seqs_to_bam_name[l["Yale-ID"]] = f'{l["Yale-ID"]}.{l["NGS_Serotype"]}.sort.bam'
                
        all_labs = set(list(seqs_to_lab.values()))


    ## prep folders
    full_genbank_paths = {}
    for i in all_labs:
        if i not in os.listdir(genbank_path):
            os.mkdir(os.path.join(genbank_path,i))

        if date_submission not in os.listdir(os.path.join(genbank_path,i)):
            os.mkdir(os.path.join(genbank_path,i,date_submission))
            
        full_genbank_paths[i] = os.path.join(genbank_path,i,date_submission)


    #put it all where it belongs

    for run, sequences in tqdm.tqdm(run_to_seqs.items()):
        
        full_file_path = os.path.join(file_path, run, "BAM")

        for seq in sequences:
            old_bam_file = os.path.join(full_file_path, seqs_to_bam_name[seq])
            
            new_bam_name = ".".join([seqs_to_bam_name[seq].split(".")[0], seqs_to_bam_name[seq].split(".")[-1]])
            new_bam_file = os.path.join(full_genbank_paths[seqs_to_lab[seq]], new_bam_name) 
            

            try:
                shutil.copy(old_bam_file, new_bam_file)
            except FileNotFoundError:
                print(old_bam_file)
                
    
if __name__ == '__main__':
    main()