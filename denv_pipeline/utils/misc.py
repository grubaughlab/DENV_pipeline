import os
import glob
import shutil

def header(v):

    print(f"""
    
            DENVseq pipeline
                Version {v}
        Verity Hill & Chrispin Chaguza
                Grubaugh Lab
    
    """)

def add_arg_to_config(key,arg,config):
    if arg:
        config[key] = arg

    return config

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