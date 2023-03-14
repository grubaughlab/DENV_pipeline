#!/usr/bin/env python3

from denv_pipeline import __version__

import os
import sys
import argparse
import snakemake
import pkg_resources

from denv_pipeline.utils import misc

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False, description=misc.header(__version__))

    parser.add_argument("--symlink", dest="symlink", help="argument for generating symlinks"),
    parser.add_argument("--run", help="number run to make folder"),

    parser.add_argument("--verbose", "-v", dest="verbose", action="store_true")
    parser.add_argument("--temp", dest="temp", help="output temporary files", action="store_true")    
    parser.add_argument("--temp-dir", dest="tempdir", help="where the temporary files go", default="temporary_files")
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
    config['temp'] = args.temp
    config['tempdir'] = args.tempdir
    config['input_file'] = args.input_file
    config["symlink"] = args.symlink
    config["cwd"] = os.path.join(cwd, args.run)
    config["denv_primers"] = pkg_resources.resource_filename('denv_pipeline', 'primers/')

    

    ## check for relevant installed stuff
    ## check for input files
    ## at the end, test if every sample ID has an associated bam file - about half way through, bam files sometimes aren't found but the script doesn't break

    # snakefile = os.path.join(thisdir,"scripts", "denv_pipeline.smk")
    # if config['verbose']:
    #     print("\n**** CONFIG ****")
    #     for k in sorted(config):
    #         print((f" - {k}: ") + f"{config[k]}")
    #     status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
    #                                 workdir=config["tempdir"],config=config,lock=False
    #                                 )
    # else:
    #     status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
    #                                 workdir=config["tempdir"],config=config,lock=False
    #                                 )

    # if status: # translate "success" into shell exit code of 0
    #     return 0

    # return 1



if __name__ == '__main__':
    main()                 
                                
                           