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

def temp_files(config, temp_files, dest):

    contains_depth = ["cons.qual.txt"]
    depth = config["depth"]
    
    for file_pattern in temp_files:
        end_pattern = ".".join(file_pattern.split(".")[1:])
        for sample in config["sample_list"]:
            for option in config["option_list"]:
                if "serotype.calls" in file_pattern:
                    if "tmp" in file_pattern:
                        name = f"tmp.{sample}.serotype.calls.{depth}.txt"
                    else:
                        name = f"{sample}.serotype.calls.txt"
                elif ".bam" in end_pattern:
                    if ".sort" in end_pattern and ".bai" not in end_pattern:
                        pass
                    else:
                        name = f"{sample}.{option}.{end_pattern}"
                else:
                    if end_pattern in contains_depth:
                        name = f"{sample}.{option}.{depth}.{end_pattern}"
                    else:
                        name = f"{sample}.{option}.{end_pattern}"

                    if config["temp"]:
                        shutil.move(os.path.join(config["outdir"], name), dest)
                    else:
                        os.remove(os.path.join(config["outdir"], name))
                

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