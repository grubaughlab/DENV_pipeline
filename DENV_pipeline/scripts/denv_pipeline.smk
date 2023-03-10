import os
import sys


rule set_up:

    input:

    output:

    parameters:

    run:
    #symlinks
    #make dirs - probably just call from misc.py
    #rename samples - also misc.py
    #ls | grep -v samples > samples.txt


rule analysis_pipeline:

    input:

    output:

    parameters:


    run:
    #denv_mapper.sh script - pull these parts out into this smk 



rule summarise_results:

    input:

    output:

    parameters:


    run:
    #DENV_summarise.sh


rule copy_to_results:

    input:

    output:

    parameters:

    run:
    #copy FINAL to  /gpfs/ycga/project/grubaugh/shared/DENVSEQ/CLINICAL/
    #make a downloadable folder (ie without BAM files) for dropbox download
