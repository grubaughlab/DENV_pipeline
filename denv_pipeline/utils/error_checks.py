import os
import sys
import csv
import pkgutil

from denv_pipeline.utils.misc import *



def check_input_files(config):

    found_sample = False

    if not os.path.exists(config["indir"]):
        sys.stderr.write(green(f"Error: cannot find samples at {config['indir']}. Please provide the input directory using `--indir`, or put them in the directory specific by `--outdir`, or (if using Yale HPC) provide the second half of the symlink provided by Yale genomics\n"))
        sys.exit(-1)
    else:
        for item in os.listdir(config["indir"]):
            if os.isdir(item):
                if os.path.join(config['indir'], item) != config["tempdir"] and item != "results" and item != "temporary_files" and item != "log_files":
                    found_sample = True
                    fastq_found = 0
                    for possible in os.listdir(os.path.join(config['indir'], item)):
                        if "fastq" in possible:
                            fastq_found += 1
                        
                            if "R1" not in possible and "R2" not in possible:
                                sys.stderr.write(green(f"Error: fastq file not named correctly. Must include 'R1' or 'R2' somewhere in the name to denote forward and reverse.\n"))
                                sys.exit(-1)
                    if fastq_found != 2:
                        sys.stderr.write(green(f"Error: Not found correct number of fastq files. There should be two in each sample folder\n"))
                        sys.exit(-1)
                
        if not found_sample:
            sys.stderr.write(green(f"Error: cannot find any valid samples at {config['indir']}. Please ensure fastq files are in directories named by sample.\n"))
            sys.exit(-1)

    return


def check_primer_dir(config):

    all_files = []
    for file in os.listdir(config["primer_directory"]):
        all_files.append(file)

    if "refs.txt" not in all_files:
        sys.stderr.write(green(f"Error: Please provide a file called 'refs.txt' containing the name of all the virus types you have references/bed files for.\n"))
        sys.exit(-1)
    
    bed = False
    trim_bed = False
    fasta_stuff = 0
    for file in all_files:
        if file.endswith(".bed"):
            if file.endswith(".trim.bed"):
                trim_bed = True
            else:
                bed = True

        if "fasta" in file:
            fasta_stuff += 1

    if not trim_bed or not bed:
        sys.stderr.write(green(f"Error: Missing bed files at {config['primer_directory']}.\n"))
        sys.exit(-1)
    if fasta_stuff != 7:
        sys.stderr.write(green(f"Error: Some reference file information missing. There should be a *.fasta, *.fasta.amb, *.fasta.ann, *.fasta.bwt, *.fasta.fai, *.fasta.pac and *.fasta.sa\n"))
        sys.exit(-1)

    return

def check_env_activated():

    if pkgutil.find_loader('snakemake') is None:
        sys.stderr.write(green(f"Error: installation not correct. Ensure that environment is activated and you have run 'pip install .'"))
        sys.exit(-1)
    

def check_ct_file(config):

    if config["ct_file"] and not config["ct_column"]:
        sys.stderr.write(green(f"Error: ct_file specified but no ct_column for ct vs coverage plot. Please provide Ct column name and id column name."))
        sys.exit(-1)

    if config["ct_file"] and not config["id_column"]:
        sys.stderr.write(green(f"Error: ct_file specified but no id_column for ct vs coverage plot. Please provide Ct column name and id column name."))
        sys.exit(-1)

    if (config["ct_column"] or config["id_column"]) and not config["ct_file"]:
        sys.stderr.write(green(f"Error: ct_column or id_column specified but no ct_file for ct vs coverage plot. Please provide file containing Ct information."))
        sys.exit(-1)

    with open(config["ct_file"]) as f:
        data = csv.DictReader(f)
        headers = data.fieldnames
        if config["ct_column"] not in headers:
            sys.stderr.write(green(f"Error: ct_column not found in ct_file"))
            sys.exit(-1)
        if config["id_column"] not in headers:
            sys.stderr.write(green(f"Error: id_column not found in ct_file"))
            sys.exit(-1)