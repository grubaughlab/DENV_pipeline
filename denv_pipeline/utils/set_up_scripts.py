import os 
import sys
import shutil
import pkg_resources


from denv_pipeline.utils.misc import *
from denv_pipeline.utils import error_checks

def overwrite(config):

    shutil.rmtree(os.path.join(config["outdir"], "results"), ignore_errors=True)
    if os.path.exists(config["tempdir"]):
        shutil.rmtree(config["tempdir"], ignore_errors=True)
    if os.path.exists(os.path.join(config["outdir"], "downloads")):
        shutil.rmtree(os.path.join(config["outdir"], "downloads"), ignore_errors=True)
    if os.path.exists(os.path.join(config["outdir"], "log_files")):
        shutil.rmtree(os.path.join(config["outdir"], "log_files"), ignore_errors=True)

def symlink_setup(config, cwd):

    symlink = config["symlink"]
    
    os.chdir(config["indir"])
    os.system(f"/home/bioinfo/software/knightlab/bin_Mar2018/ycgaFastq {symlink}")
    os.chdir(cwd)
    
    for sample_dir in os.listdir(config["indir"]):
        upper_path = os.path.join(config["indir"], sample_dir)
        if os.path.isdir(upper_path) and sample_dir != "results" and os.path.join(config["indir"], sample_dir) != config["tempdir"]:
            lower_path = os.path.join(upper_path, "Unaligned")
            for file in os.listdir(lower_path):
                if file.endswith(".fastq.gz"):
                    shutil.move(os.path.join(lower_path, file), upper_path)

            shutil.rmtree(lower_path)

    return config


def get_sample_list(config):

    config["sample_list"] = []
    for sample_dir in os.listdir(config["indir"]):
        if os.path.isdir(os.path.join(config["indir"], sample_dir)):
            if sample_dir != "results" and os.path.join(config["indir"], sample_dir) != config["tempdir"] and sample_dir != "log_files":
                config["sample_list"].append(sample_dir)

    return config


def make_folders(config):

    outdir = config["outdir"]

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    if os.path.exists(os.path.join(outdir, "results")) and not config["overwrite"]:
        sys.stderr.write(green(f"Error: results files already exist at {outdir}. Use --overwrite flag to delete and regenerate results."))
        sys.exit(-1)
    else:
        os.mkdir(os.path.join(outdir, "results"))
        if config["temp"]:
            if not os.path.exists(config["tempdir"]):
                os.mkdir(config["tempdir"])

    if not os.path.exists(os.path.join(config["outdir"], "log_files")):
        os.mkdir(os.path.join(config["outdir"], "log_files"))

def set_up_primer_directory(config, args):

    if args.primer_directory:
        if not args.primer_directory.endswith("/"):
            config["primer_directory"] = f'{args.primer_directory}/'
        else:
            config["primer_directory"] = args.primer_directory
        error_checks.check_primer_dir(config)
    else:
        if config["verbose"]:
            print("Using DENV primers")
        config["primer_directory"] = pkg_resources.resource_filename('denv_pipeline', 'DENV_primers_and_refs/')

    config["virus_type_list"] = []
    with open(os.path.join(config["primer_directory"], "refs.txt")) as f:
        for l in f:
            config["virus_type_list"].append(l.strip("\n"))

    return config
    


def find_fastq_names(config):

    config["fastq_R1"] = {}
    config["fastq_R2"] = {}
    for sample in config["sample_list"]:
        for file in os.listdir(os.path.join(config["indir"], sample)):
           if "fastq" in file and "R1" in file:
               config["fastq_R1"][sample] = file
               config["fastq_R2"][sample] = file.replace("R1", "R2")

    return config
                   
            
