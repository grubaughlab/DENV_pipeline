#!/usr/bin/env python3

from denv_pipeline import __version__

import os
import sys
import argparse
import snakemake
import pkg_resources
import datetime as dt

from denv_pipeline.utils import misc

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False, description=misc.header(__version__))

    parser.add_argument("--symlink", dest="symlink", help="argument for generating symlinks", default="")
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
    
    if not args.outdir:
        outdir = f'denv_seq_{dt.datetime.today().date()}'
    else:
        outdir = (args.outdir).rstrip("/")
    
    config["outdir"] = outdir
    config["tempdir"] = os.path.join(outdir, args.tempdir)

    if not args.indir:
        config["indir"] = config["outdir"]
    else:
        config["indir"] = args.indir


    if args.primer_directory:
        config["primer_directory"] = args.primer_directory
    else:
        if config["verbose"]:
            print("Using DENV primers")
            config["primer_directory"] = pkg_resources.resource_filename('denv_pipeline', 'primers/')
    

    misc.check_input_files(config)

    config["sample_list"] = []
    for sample_dir in os.listdir(config["indir"]):
        if os.path.isdir(os.path.join(config["indir"], sample_dir)):
            if sample_dir != "download" and os.path.join(config["indir"], sample_dir) != config["tempdir"]:
                config["sample_list"].append(sample_dir)

    config["option_list"] = []
    with open(os.path.join(config["primer_directory"], "refs.txt")) as f:
        for l in f:
            config["option_list"].append(l.strip("\n"))

    if config["overwrite"]:
        misc.overwrite(config)
        
    
    ## check for relevant installed stuff
    ## check for input files - either the symlink is present, or if an "indir" is used then they should be in there already. Also in the right format
    

    snakefile = os.path.join(thisdir,"scripts", "denv_pipeline.smk")
    if config['verbose']:
        print("\n**** CONFIG ****")
        for k in sorted(config):
            print((f" - {k}: ") + f"{config[k]}")
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True, force_incomplete=True,
                                    workdir=cwd,config=config,lock=False
                                    )
    else:
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True, force_incomplete=True,
                                    workdir=cwd,config=config,lock=False
                                    )
        

    #QC plots

    if status: # translate "success" into shell exit code of 0
        return 0

    return 1



if __name__ == '__main__':
    main()                 
                                
                           