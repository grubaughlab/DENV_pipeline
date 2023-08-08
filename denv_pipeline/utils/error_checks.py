import os
import sys
import csv
import pkgutil
import re
import yaml

from denv_pipeline.utils.misc import *
from denv_pipeline.utils import set_up_scripts

def check_configfile(cwd,config_arg):
    
    configfile = os.path.join(cwd,config_arg)
    ending = configfile.split(".")[-1]

    if ending not in ["yaml","yml"]:
        sys.stderr.write(cyan(f'Error: config file {configfile} must be in yaml format.\n'))
        sys.exit(-1)
    
    elif not os.path.isfile(configfile):
        sys.stderr.write(cyan(f'Error: cannot find config file at {configfile}\n'))
        sys.exit(-1)
   
    else:    
        print(green(f"Input config file:") + f" {configfile}")
        return configfile
    
def fix_paths(config):

    path_args = ["outdir", "indir", "tempdir", "reference_directory", "config"]
    
    for arg in path_args:
        if config[arg] and " " in config[arg]:
            sys.stderr.write(green(f"Error: space in {arg} path name. Please change directory to where there is no spaces. NB dropbox paths always have a spcae in\n"))
            sys.exit(-1)
        
    return config


def check_input_files(config):

    found_sample = False

    if not os.path.exists(config['indir']):
        sys.stderr.write(green(f"Error: cannot find {config['indir']}. Please provide the input directory using `--indir`, or put them in the directory specific by `--outdir`, or (if using Yale HPC) provide the second half of the symlink provided by Yale genomics\n"))
        sys.exit(-1)
    else:
        for item in os.listdir(config["indir"]):
            if os.path.isdir(os.path.join(config["indir"], item)):
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

    if not os.path.exists(config["reference_directory"]):
        sys.stderr.write(green(f"Error: reference directory not found at {config['reference_directory']}"))
        sys.exit(-1)

    all_files = []
    for file in os.listdir(config["reference_directory"]):
        all_files.append(file)

    if "refs.txt" not in all_files:
        sys.stderr.write(green(f"Error: Please provide a file called 'refs.txt' containing the name of all the virus types you have references/bed files for.\n"))
        sys.exit(-1)
    else:
        virus_types = []
        with open(os.path.join(config["reference_directory"], "refs.txt")) as f:
            for l in f:
                virus_types.append(l.strip("\n"))

        for virus_type in virus_types:
            if virus_type != "":

                bed_file = f'{virus_type}.bed'
                fasta_file = f'{virus_type}.fasta'

                if not bed_file in all_files:
                    sys.stderr.write(green(f"Error: Missing bed file for {virus_type} in {config['reference_directory']}. I am expecting it to be called {bed_file} and the '.bed' has to be lower case.\n"))
                    sys.exit(-1)
                if not fasta_file in all_files:
                    sys.stderr.write(green(f"Error: Missing reference file for {virus_type} in {config['reference_directory']}. I am expecting it to be called {fasta_file}\n"))
                    sys.exit(-1)
            
                with open(os.path.join(config["reference_directory"],fasta_file)) as f:
                    seq_count = 0
                    for l in f:
                        if l.startswith(">"):
                            seq_count += 1
                if seq_count != 1:
                    sys.stderr.write(green(f"Error: Wrong number of sequences in reference file {fasta_file} - there should only be one.\n"))
                    sys.exit(-1)
            
            else:
                sys.stderr.write(green(f"Error: Empty line in refs.txt. Please remove.\n"))
                sys.exit(-1)

        return

def check_env_activated():

    if pkgutil.find_loader('snakemake') is None or pkgutil.find_loader('Bio') is None:
        sys.stderr.write(green(f"Error: installation not correct. Ensure that environment is activated and you have run 'pip install .'\n"))
        sys.exit(-1)
    

def check_ct_file(config):

    if config["ct_file"] and not config["ct_column"]:
        sys.stderr.write(green(f"Error: ct_file specified but no ct_column for ct vs coverage plot. Please provide Ct column name and id column name.\n"))
        sys.exit(-1)

    if config["ct_file"] and not config["id_column"]:
        sys.stderr.write(green(f"Error: ct_file specified but no id_column for ct vs coverage plot. Please provide Ct column name and id column name.\n"))
        sys.exit(-1)

    if (config["ct_column"] or config["id_column"]) and not config["ct_file"]:
        sys.stderr.write(green(f"Error: ct_column or id_column specified but no ct_file for ct vs coverage plot. Please provide file containing Ct information.\n"))
        sys.exit(-1)

    if not os.path.exists(config["ct_file"]):
        sys.stderr.write(green(f"Error: Ct file not found at {config['ct_file']}"))
        sys.exit(-1)


    with open(config["ct_file"]) as f:
        data = csv.DictReader(f)
        headers = data.fieldnames
        if config["ct_column"] not in headers:
            sys.stderr.write(green(f"Error: {config['ct_column']} not found in ct_file\n"))
            sys.exit(-1)
        if config["id_column"] not in headers:
            sys.stderr.write(green(f"Error: {config['id_column']} not found in ct_file\n"))
            sys.exit(-1)


def check_threshold(config):

    if float(config["threshold"]) > 1:
        sys.stderr.write(green(f"Error: consensus threshold must be between 0 and 1\n"))
        sys.exit(-1)