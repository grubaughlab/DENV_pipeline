import os
import sys
import csv

from denv_pipeline.utils.misc import *



def check_input_files(config):

    #if input != output, does it exist (can put symlinks in diff place so not reliant on that)
    #if not symlink, does it exist and is the structure correct - ie sample_name folder with 2xfastq files in it with R1 and R2 on them

    return


def check_ct_file(config):

    if config["ct_file"] and not config["ct_column"]:
        sys.stderr.write(green(f"Error: ct_file specified but no ct_column for ct vs coverage plot. Please provide Ct column name and id column name."))
        sys.exit(-1)

    if config["ct_file"] and not config["id_column"]:
        sys.stderr.write(green(f"Error: ct_file specified but no id_column for ct vs coverage plot. Please provide Ct column name and id column name."))
        sys.exit(-1)

    if (config["ct_column"] or config["id_column"]) and not config["ct_file"]:
        sys.stderr.write(green(f"Error: ct_column or id_column specified but no ct_file for ct vs coverage plot. Please provide file containing Ct information."))
        sys.exit(-1)

    with open(config["ct_file"]) as f:
        data = csv.DictReader(f)
        headers = data.fieldnames
        if config["ct_column"] not in headers:
            sys.stderr.write(green(f"Error: ct_column not found in ct_file"))
            sys.exit(-1)
        if config["id_column"] not in headers:
            sys.stderr.write(green(f"Error: id_column not found in ct_file"))
            sys.exit(-1)