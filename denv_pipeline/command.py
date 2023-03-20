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
    
    parser.add_argument("--temp", dest="temp", action="store_true", help="keep intermediate files")
    parser.add_argument("--tempdir", dest="tempdir", help="where the temporary files go", default="temporary_files")
    parser.add_argument("--download", action="store_true", help="make a folder without bam files for download")

    parser.add_argument("--slurm", help="flag for if running on HPC with slurm", action="store_true")
    parser.add_argument("--verbose", "-v", dest="verbose", action="store_true")
    parser.add_argument("--help", "-h", action="store_true", dest="help")


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
    config["indir"] = args.indir
    config["slurm"] = args.slurm
    config["temp"] = args.temp
    config["download"] = args.download
    
    config["denv_primers"] = pkg_resources.resource_filename('denv_pipeline', 'primers/')

    if not args.outdir:
        outdir = f'denv_seq_{dt.datetime.today().date()}'
    else:
        outdir = (args.outdir).rstrip("/")
    
    config["outdir"] = outdir
    config["tempdir"] = os.path.join(outdir, args.tempdir)

    if not args.indir:
        config["indir"] = config["outdir"]
    

    misc.check_input_files(config)
    
    ## check for relevant installed stuff
    ## check for input files - either the symlink is present, or if an "indir" is used then they should be in there already. Also in the right format
    

    snakefile = os.path.join(thisdir,"scripts", "denv_pipeline.smk")
    if config['verbose']:
        print("\n**** CONFIG ****")
        for k in sorted(config):
            print((f" - {k}: ") + f"{config[k]}")
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=config["cwd"],config=config,lock=False
                                    )
    else:
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=config["cwd"],config=config,lock=False
                                    )
        


    #denv_summarise
    #QC plots
    #copy to right place - if temp tru, put in temporary folder, if not then delete it all

    if status: # translate "success" into shell exit code of 0
        return 0

    return 1



if __name__ == '__main__':
    main()                 
                                
                           