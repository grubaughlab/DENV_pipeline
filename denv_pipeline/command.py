#!/usr/bin/env python3

from denv_pipeline import __version__

from denv_pipeline.utils import misc
from denv_pipeline.utils import error_checks
from denv_pipeline.utils import set_up_scripts
import pkg_resources
import os
import sys
import argparse
import datetime as dt
import yaml

error_checks.check_env_activated()

import snakemake


cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False, description=misc.header(__version__))

    parser.add_argument("--config", help="config file containing all relevant arguments")
    parser.add_argument("--dry-run", dest="dry_run", action="store_true", help="do all error checks and make files but don't run the pipeline.")

    parser.add_argument("--symlink", dest="symlink", help="argument for generating symlinks")
    parser.add_argument("--indir", help="directory containing samples. Each sample must be a folder with the forward and reverse runs in. Default is same as output directory")
    parser.add_argument("--outdir", dest="outdir", help="location where files will be stored.")
    parser.add_argument("--reference-directory", "-rd", help="location where bed files and reference genomes are")
    parser.add_argument("--depth", help="depth to map sequences to. Default=10")
    parser.add_argument("--threshold", help="threshold to call consensus positions at, default=0.75",dest="threshold")
    parser.add_argument("--no-cap", help="remove cap from number of reads used for consensus generation. Default uses 10,000 reads.",dest="no_cap", action="store_true")
    
    parser.add_argument("--temp", dest="temp", action="store_true", help="keep intermediate files")
    parser.add_argument("--tempdir", dest="tempdir", help="where the temporary files go")
    parser.add_argument("--download", action="store_true", help="make a folder without bam files for download")

    parser.add_argument("--slurm", help="flag for if running on HPC with slurm", action="store_true")
    parser.add_argument("--slurm-cores", help="number of slurm cores to assign. Default is 10", dest="slurm_cores", type=int)
    parser.add_argument("--cores", help="number of non-slurm cores to assign. Default is 1", type=int)
    parser.add_argument("--verbose", "-v", dest="verbose", action="store_true")
    parser.add_argument("--help", "-h", action="store_true", dest="help")
    parser.add_argument("--overwrite", help="overwrite current results", action="store_true")

    parser.add_argument("--ct-file",dest="ct_file", help="to produce a plot of Ct against coverage, provide a csv file containing Ct information by sample")
    parser.add_argument("--ct-column", dest="ct_column", help="Name of Ct column in Ct file for plot")
    parser.add_argument("--id-column", dest="id_column", help="Name of ID column in Ct file to make Ct plot")

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)


    config = {}
    config = set_up_scripts.get_defaults(config)

    if args.config:
        configfile = error_checks.check_configfile(cwd,args.config)
        config = set_up_scripts.parse_yaml_file(configfile,config)
        config["config"] = True
    else:
        configfile = False

    for arg_name, arg_value in vars(args).items():
        config = misc.add_arg_to_config(arg_name, arg_value, config)
    if "indir" not in config:
        config["indir"] = config["outdir"]

    error_checks.check_threshold(config)

    config["outdir"] = config["outdir"].rstrip("/")
    config = set_up_scripts.set_up_temporary_directory_path(config)
    config = set_up_scripts.set_up_reference_directory(config)

    if config["overwrite"]:
        set_up_scripts.overwrite(config)
    set_up_scripts.make_folders(config)

    if config["symlink"]:
        config = set_up_scripts.symlink_setup(config, cwd)
    else:
        error_checks.check_input_files(config)

    config = set_up_scripts.get_sample_list(config)
    
    if config["ct_file"] or config["ct_column"] or config["id_column"]:
        error_checks.check_ct_file(config)    

    if config["no_cap"]:
        config["cap"] = 0
    else:
        config["cap"] = 10000

    set_up_scripts.output_config(config)

    snakefile = os.path.join(thisdir,"scripts", "denv_pipeline.smk")
    if config['verbose'] or config["dry_run"]:
        print("\n**** CONFIG ****")
        for k in sorted(config):
            print((f" - {k}: ") + f"{config[k]}")

        num_samples = len(config["sample_list"])
        print(f"\n Analysing {num_samples} samples against reference files for: \n")
        for k in sorted(config["virus_type_list"]):
            print(f" - {k} ")
   
    if not config["dry_run"]:
        if config["slurm"]:
            status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True, force_incomplete=True, latency_wait=10,
                                    workdir=cwd,config=config,lock=False, slurm=True, cores=config["slurm_cores"]
                                    )
        elif config["cores"]:
            status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True, force_incomplete=True,
                                    workdir=cwd,config=config,lock=False, cores=config["cores"]
            )
        else:
            status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=cwd,config=config,lock=False
                                    )


        if status: # translate "success" into shell exit code of 0
            return 0

        return 1

    else:
        sys.stderr.write("Dry run complete")
        return 0


if __name__ == '__main__':
    main()                 
                                
                           
