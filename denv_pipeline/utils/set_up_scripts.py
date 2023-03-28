import os 
import sys
import shutil


from denv_pipeline.utils.misc import *


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
            if os.path.join(config["indir"], sample_dir) != "results" and os.path.join(config["indir"], sample_dir) != config["tempdir"]:
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
