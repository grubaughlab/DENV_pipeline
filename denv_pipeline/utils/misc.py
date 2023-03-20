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


def check_input_files(config):

    #if input != output, does it exist (can put symlinks in diff place so not reliant on that)
    #if not symlink, does it exist and is the structure correct - ie sample_name folder with 2xfastq files in it with R1 and R2 on them

    return