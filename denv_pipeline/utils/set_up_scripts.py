import os 
import sys
import shutil
import pkg_resources
import yaml
import datetime as dt

from denv_pipeline.utils import misc
from denv_pipeline.utils import error_checks

def get_defaults(config):
       
    config["config"] = False
    
    config["depth"] = 10
    config["threshold"] = 0.75
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
    config["dry_run"] = False

    config["outdir"] = f'seq_analysis_{dt.datetime.today().date()}' 

    config["reference_directory"] = pkg_resources.resource_filename('denv_pipeline', 'DENV_primers_and_refs')

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
        sys.stderr.write(misc.cyan(f'Error: invalid key in config file.\n') + f'\t- {invalid_keys[0]}\n')
        sys.exit(-1)
    elif len(invalid_keys) >1:
        keys = ""
        for i in invalid_keys:
            keys += f"\t- {i}\n"
        sys.stderr.write(misc.cyan(f'Error: invalid keys in config file.\n') + f'{keys}')
        sys.exit(-1)
    print(misc.green(f"Adding {overwriting} arguments to internal config."))

    return configdict

def load_yaml(f):
    try:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
    except:
        sys.stderr.write(misc.cyan(f'Error: failed to read config file. Ensure your file in correct yaml format.\n'))
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

    non_sample_dirs = ["results", "log_files", "downloads"]

    symlink = config["symlink"]
    
    os.chdir(config["indir"])
    os.system(f"ycgaFastq {symlink}")
    os.chdir(cwd)
    
    for sample_dir in os.listdir(config["indir"]):
        upper_path = os.path.join(config["indir"], sample_dir)
        if os.path.isdir(upper_path) and sample_dir not in non_sample_dirs and os.path.join(config["indir"], sample_dir) != config["tempdir"]:
            lower_path = os.path.join(upper_path, "Unaligned")
            for file in os.listdir(lower_path):
                if file.endswith(".fastq.gz"):
                    if not os.path.exists(os.path.join(upper_path, file)):
                        shutil.move(os.path.join(lower_path, file), upper_path)

            shutil.rmtree(lower_path)

    return config


def get_sample_list(config):

    non_sample_dirs = ["results", "log_files", "downloads"]

    config["sample_list"] = []
    for sample_dir in os.listdir(config["indir"]):
        if os.path.isdir(os.path.join(config["indir"], sample_dir)):
            if sample_dir not in non_sample_dirs and os.path.join(config["indir"], sample_dir) != config["tempdir"]:
                config["sample_list"].append(sample_dir)

    return config


def make_folders(config):

    out_dir = config["outdir"]
    temp_dir = config["tempdir"]
    results_dir = os.path.join(config["outdir"], "results")
    log_dir = os.path.join(config["outdir"], "log_files")

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    if os.path.exists(os.path.join(out_dir, "results")) and not config["overwrite"]:
        sys.stderr.write(misc.green(f"Error: results files already exist at {out_dir}. Use --overwrite flag to delete and regenerate results.\n"))
        sys.exit(-1)
    else:
        misc.make_directory(results_dir)
        misc.make_directory(temp_dir)
        misc.make_directory(log_dir)

        misc.make_directory(os.path.join(results_dir, "bam_files"))
        misc.make_directory(os.path.join(results_dir, "variants"))
        misc.make_directory(os.path.join(results_dir, "depth"))
        misc.make_directory(os.path.join(results_dir, "consensus_sequences"))
        misc.make_directory(os.path.join(results_dir, "alignments"))

        if config["download"]:
            misc.make_directory(os.path.join(config["outdir"], "downloads")) 

def set_up_reference_directory(config):

    reference_directory = config["reference_directory"]

    if reference_directory != pkg_resources.resource_filename('denv_pipeline', 'DENV_primers_and_refs'):
        config['reference_directory'] = reference_directory.rstrip("/")
        error_checks.check_primer_dir(config)
        
    config["virus_type_list"] = []
    with open(os.path.join(config["reference_directory"], "refs.txt")) as f:
        for l in f:
            config["virus_type_list"].append(l.strip("\n"))

    return config

def set_up_temporary_directory_path(config):

    if config["outdir"] not in config["tempdir"]:
        temp = os.path.join(config["outdir"], config["tempdir"])
    else:
        temp = config["tempdir"]

    config["tempdir"] = temp.rstrip("/")

    return config

    
       
def get_valid_keys():

    valid_keys = []

    valid_keys.append("verbose")
    valid_keys.append("symlink")
    valid_keys.append("indir")
    valid_keys.append("outdir")
    valid_keys.append("reference_directory")
    valid_keys.append("depth")
    valid_keys.append("threshold")
    valid_keys.append("temp")
    valid_keys.append("tempdir")
    valid_keys.append("download")
    valid_keys.append("slurm")
    valid_keys.append("slurm_cores")
    valid_keys.append("verbose")
    valid_keys.append("config")
    valid_keys.append("dry_run")
    valid_keys.append("overwrite")
    valid_keys.append("ct_file")
    valid_keys.append("ct_column")
    valid_keys.append("id_column")

    return valid_keys

def output_config(config):

    valid_keys = get_valid_keys()

    new_config = {}
    for k,v in config.items():
        if k in valid_keys:
            new_config[k] = v

    with open(os.path.join(config["outdir"], "output_config.yml"), 'w') as file:
        documents = yaml.dump(new_config, file)