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
        os.remove(file)

def remove_multiple_files(pattern):
    for f in glob.glob(pattern):
        os.remove(f)