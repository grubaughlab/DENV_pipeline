import os
import glob
import shutil

def header(v):

    print(f"""
    
            DENVseq pipeline
                Version {v}
        Chrispin Chaguza & Verity Hill
                Grubaugh Lab
    
    """)


def remove_file(file):
    if os.path.exists(file):
        print(f"removing {file}")
        os.remove(file)

def remove_multiple_files(pattern):
    for f in glob.glob(pattern):
        print(f"removing {pattern}")
        os.remove(f)

def make_directory(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)


def overwrite(config):

    shutil.rmtree(os.path.join(config["outdir"], "results"), ignore_errors=True)
    if os.path.exists(config["tempdir"]):
        shutil.rmtree(config["tempdir"], ignore_errors=True)
    if os.path.exists(os.path.join(config["outdir"], "downloads")):
        shutil.rmtree(os.path.join(config["outdir"], "downloads"), ignore_errors=True)
    
    if os.path.exists(os.path.join(config["outdir"], "samples.txt")):
        os.remove(os.path.join(config["outdir"], "samples.txt"))
    if os.path.exists(os.path.join(config["outdir"], "jobs.txt")):
        os.remove(os.path.join(config["outdir"], "jobs.txt"))


def temp_files(config, temp_files, dest):

    contains_depth = ["cons.qual.txt", "variants.tsv", "out.aln", "out.trim.aln"]
    depth = config["depth"]
    for file_pattern in temp_files:
        for sample in config["sample_list"]:
            if "serotype.calls" in file_pattern:
                if "tmp" in file_pattern:
                    name = f"tmp.{sample}.serotype.calls.{depth}.txt"
                else:
                    name = f"{sample}.serotype.calls.txt"

                if config["temp"]:
                    shutil.move(os.path.join(config["outdir"], name), dest)
                else:
                    os.remove(os.path.join(config["outdir"], name))
            
            else:
                for option in config["option_list"]:
                    if ".bam" in file_pattern:
                        if ".sort" in file_pattern and ".bai" not in file_pattern:
                            pass
                        else:
                            name = f"{sample}.{option}.{file_pattern}"
                    else:
                        if file_pattern in contains_depth:
                            name = f"{sample}.{option}.{depth}.{file_pattern}"
                        else:
                            name = f"{sample}.{option}.{file_pattern}"

                    if config["temp"]:
                        shutil.move(os.path.join(config["outdir"], name), dest)
                    else:
                        os.remove(os.path.join(config["outdir"], name))

    

def alignment_components(config, temp_dir):

    path = os.path.join(config["outdir"],"results", "alignment")
    depth = config["depth"]
    components = []
    for option in config["option_list"]:
        for sample in config["sample_list"]:
            file1 = f'{sample}.{option}.{depth}.out.aln'
            file2 = f'{sample}.{option}.{depth}.out.trim.aln'

            if os.path.exists(os.path.join(path, file1)):
                if config["temp"]:
                    shutil.move(os.path.join(path, file1), temp_dir)
                    shutil.move(os.path.join(path, file2), temp_dir)
                else:
                    os.remove(os.path.join(path, file1))
                    os.remove(os.path.join(path, file2))
                

def check_input_files(config):

    #if input != output, does it exist (can put symlinks in diff place so not reliant on that)
    #if not symlink, does it exist and is the structure correct - ie sample_name folder with 2xfastq files in it with R1 and R2 on them

    return

END_FORMATTING = '\033[0m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'

def red(text):
    return RED + text + END_FORMATTING

def cyan(text):
    return CYAN + text + END_FORMATTING

def green(text):
    return GREEN + text + END_FORMATTING

def yellow(text):
    return YELLOW + text + END_FORMATTING