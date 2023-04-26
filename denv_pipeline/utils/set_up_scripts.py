import os 
import sys
import shutil
import pkg_resources
import yaml
import datetime as dt

from denv_pipeline.utils.misc import *
from denv_pipeline.utils import error_checks

def get_defaults(config):

    config["depth"] = 20
    config["tempdir"] = "temporary_files"
    
    config['ct_file'] = False
    config["ct_column"] = False
    config["id_column"] = False
    
    config["verbose"] = False
    config["slurm"] = False
    config["slurm_cores"] = 10
    config["download"] = False
    config["temp"] = False
    config["overwrite"] = False
    config["symlink"] = False

    config["outdir"] = f'seq_analysis_{dt.datetime.today().date()}' 
    config["indir"] = config["outdir"]

    config["reference_directory"] = pkg_resources.resource_filename('denv_pipeline', 'DENV_primers_and_refs/')

    return config

def parse_yaml_file(configfile,configdict):
    
    overwriting = 0
    path_to_file = os.path.abspath(os.path.dirname(configfile))
    invalid_keys = []
    
    valid_keys = get_valid_keys()

    with open(configfile,"r") as f:
        input_config = load_yaml(f) # try load file else exit with msg

        for key in input_config:
            value = input_config[key]
            if value == None: # dont count blank entries
                pass
            else:
                clean_key = key.lstrip("-").replace("-","_").rstrip(" ").lstrip(" ").lower()

                if not clean_key in valid_keys:
                    invalid_keys.append(key)
                    break
                    
                configdict[clean_key] = value
                overwriting += 1

    if len(invalid_keys)==1:
        sys.stderr.write(cyan(f'Error: invalid key in config file.\n') + f'\t- {invalid_keys[0]}\n')
        sys.exit(-1)
    elif len(invalid_keys) >1:
        keys = ""
        for i in invalid_keys:
            keys += f"\t- {i}\n"
        sys.stderr.write(cyan(f'Error: invalid keys in config file.\n') + f'{keys}')
        sys.exit(-1)
    print(green(f"Adding {overwriting} arguments to internal config."))

    return configdict

def load_yaml(f):
    try:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
    except:
        sys.stderr.write(cyan(f'Error: failed to read config file. Ensure your file in correct yaml format.\n'))
        sys.exit(-1)
    return input_config


def overwrite(config):

    shutil.rmtree(os.path.join(config["outdir"], "results"), ignore_errors=True)
    if os.path.exists(config["tempdir"]):
        shutil.rmtree(config["tempdir"], ignore_errors=True)
    if os.path.exists(os.path.join(config["outdir"], "downloads")):
        shutil.rmtree(os.path.join(config["outdir"], "downloads"), ignore_errors=True)
    if os.path.exists(os.path.join(config["outdir"], "log_files")):
        shutil.rmtree(os.path.join(config["outdir"], "log_files"), ignore_errors=True)
    if os.path.exists(os.path.join(config["outdir"], ".snakemake")):
        shutil.rmtree(os.path.join(config["outdir"], ".snakemake"), ignore_errors=True)
    

def symlink_setup(config, cwd):

    symlink = config["symlink"]
    
    os.chdir(config["indir"])
    os.system(f"/home/bioinfo/software/knightlab/bin_Mar2018/ycgaFastq {symlink}")
    os.chdir(cwd)
    
    for sample_dir in os.listdir(config["indir"]):
        upper_path = os.path.join(config["indir"], sample_dir)
        if os.path.isdir(upper_path) and sample_dir != "results" and sample_dir != "log_files" and os.path.join(config["indir"], sample_dir) != config["tempdir"]:
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

def set_up_reference_directory(config):

    reference_directory = config["reference_directory"]

    if reference_directory:
        if not reference_directory.endswith("/"):
            config["reference_directory"] = f'{reference_directory}/'
        error_checks.check_primer_dir(config)
        
    config["virus_type_list"] = []
    with open(os.path.join(config["reference_directory"], "refs.txt")) as f:
        for l in f:
            config["virus_type_list"].append(l.strip("\n"))

    return config
    
       
def get_valid_keys():

    valid_keys = []

    valid_keys.append("verbose")
    valid_keys.append("symlink")
    valid_keys.append("indir")
    valid_keys.append("outdir")
    valid_keys.append("reference_directory")
    valid_keys.append("depth")
    valid_keys.append("temp")
    valid_keys.append("tempdir")
    valid_keys.append("download")
    valid_keys.append("slurm")
    valid_keys.append("slurm_cores")
    valid_keys.append("verbose")
    valid_keys.append("overwrite")
    valid_keys.append("ct_file")
    valid_keys.append("ct_column")
    valid_keys.append("id_column")

    return valid_keys