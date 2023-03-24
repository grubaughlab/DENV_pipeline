import os
import glob

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