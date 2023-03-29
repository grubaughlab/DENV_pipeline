#!/usr/bin/env python3

from denv_pipeline import __version__

import os
import sys
import argparse
import snakemake
import pkg_resources
import datetime as dt

from denv_pipeline.utils import misc
from denv_pipeline.utils import error_checks
from denv_pipeline.utils import set_up_scripts

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False, description=misc.header(__version__))

    parser.add_argument("--symlink", dest="symlink", help="argument for generating symlinks", default=None)
    parser.add_argument("--indir", help="directory containing samples. Each sample must be a folder with the forward and reverse runs in. Default is same as output directory")
    parser.add_argument("--outdir", dest="outdir", help="location where files will be stored.")
    parser.add_argument("--primer-directory", "-pd", help="location where bed files etc for references are")
    parser.add_argument("--depth", help="depth to map sequences to. Default=20", default=20)
    
    parser.add_argument("--temp", dest="temp", action="store_true", help="keep intermediate files")
    parser.add_argument("--tempdir", dest="tempdir", help="where the temporary files go", default="temporary_files")
    parser.add_argument("--download", action="store_true", help="make a folder without bam files for download")

    parser.add_argument("--slurm", help="flag for if running on HPC with slurm", action="store_true")
    parser.add_argument("--verbose", "-v", dest="verbose", action="store_true")
    parser.add_argument("--help", "-h", action="store_true", dest="help")
    parser.add_argument("--overwrite", help="overwrite current results", action="store_true")

    parser.add_argument("--ct-file",dest="ct_file", help="to produce a plot of Ct against coverage, provide a csv file containing Ct information by sample", default=None)
    parser.add_argument("--ct-column", dest="ct_column", help="Name of Ct column in Ct file for plot", default=None)
    parser.add_argument("--id-column", dest="id_column", help="Name of ID column in Ct file to make Ct plot", default=None)

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)


    config = {}
    config['verbose'] = args.verbose
    config["symlink"] = args.symlink
    config["slurm"] = args.slurm
    config["temp"] = args.temp
    config["download"] = args.download
    config["overwrite"] = args.overwrite
    config["depth"] = args.depth
    config["ct_file"] = args.ct_file
    config["ct_column"] = args.ct_column
    config["id_column"] = args.id_column
    
    if not args.outdir:
        outdir = f'seq_analysis_{dt.datetime.today().date()}'
    else:
        outdir = (args.outdir).rstrip("/")
    
    config["outdir"] = outdir
    config["tempdir"] = os.path.join(outdir, args.tempdir)

    if not args.indir:
        config["indir"] = config["outdir"]
    else:
        config["indir"] = args.indir

    if config["overwrite"]:
        misc.overwrite(config)

    set_up_scripts.make_folders(config)
    
    if args.symlink:
        config = set_up_scripts.symlink_setup(config, cwd)

    config = set_up_scripts.get_sample_list(config)
    conifg = set_up_scripts.find_fastq_names(config)
    error_checks.check_input_files(config)
    
    if args.primer_directory:
        config["primer_directory"] = args.primer_directory
    else:
        if config["verbose"]:
            print("Using DENV primers")
        config["primer_directory"] = pkg_resources.resource_filename('denv_pipeline', 'primers/')

    config["option_list"] = []
    with open(os.path.join(config["primer_directory"], "refs.txt")) as f:
        for l in f:
            config["option_list"].append(l.strip("\n"))

    if config["ct_file"] or config["ct_column"] or config["id_column"]:
        error_checks.check_ct_file(config)
        
    
    ## check for relevant installed stuff
    ## check for input files - either the symlink is present, or if an "indir" is used then they should be in there already. Also in the right format
    

    snakefile = os.path.join(thisdir,"scripts", "denv_pipeline.smk")
    if config['verbose']:
        print("\n**** CONFIG ****")
        for k in sorted(config):
            print((f" - {k}: ") + f"{config[k]}")
    if config["slurm"]:
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                workdir=cwd,config=config,lock=False, slurm=True
                                )
    else:
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                workdir=cwd,config=config,lock=False
                                )


    #QC plots

    if status: # translate "success" into shell exit code of 0
        return 0

    return 1



if __name__ == '__main__':
    main()                 
                                
                           
