import os
import sys
import datetime as dt

from misc import *

rule all: #this will be the outputs we want
    input:
        os.path.join(config["cwd"], "jobs.txt")

rule setup:
    output:
        sample_file = os.path.join(config["cwd"], "samples.txt")
    params:
        cwd = config["cwd"],
        symlink = config["symlink"]
    run:
        if not os.path.exists(config["cwd"]):
            os.mkdir(config["cwd"])
        
        shell("cd {params.cwd}")
        
        #snakemake force rerun or something?
        misc.remove_file("DENV.serotype.calls.tsv")
        misc.remove_multiple_files("*.serotype.txt")
        misc.remove_multiple_files("tmp.*.serotype.calls.*.txt")
        misc.remove_multiple_files("*.*.out.aln")

        #shell("/home/bioinfo/software/knightlab/bin_Mar2018/ycgaFastq {params.symlink:q}")
        shell("ls | grep -v samples > {output.sample_file:q}")


rule denv_mapper:
    input:
        mapper_script = os.path.join(workflow.current_basedir,"DENV_MAPPER.sh"),
        sample_file = rules.setup.output.sample_file,
        refs = os.path.join(config["denv_primers"], "DENV.refs.txt")
        python_script = os.path.join(workflow.current_basedir,"serotypeCaller.py")
    output:
        jobs = os.path.join(config["cwd"], "jobs.txt")
        #have the full denv ones here
    run:
        with open(output.jobs, 'w') as fw:
            with open(input.sample_file) as f:
                for l in f:
                    name = l.strip("\n")
                    basename = name.split("_")[0]
                    fw.write(f"bash {input.mapper_script} {basename} {name}/*/{name}*_R1_*.fastq.gz {name}/*/{name}*_R2_*.fastq.gz {input.refs} {input.python_script}")

            
        if config["slurm"]:
            print("preparing for slurm run")
            shell("dsq --job-name denv.mapper --job-file {output.jobs:q} --mem-per-cpu=10G --cpus-per-task=1") 
            filename = f"dsq-jobs-{dt.datetime.today().date()}.sh"
            shell("sbatch {filename}")
        
        else:
            print("running each sample sequentially")
            with open(output.jobs) as f:
                for l in f:
                    command = l.strip("\n")
                    shell("{command}")
        
            

    rule denv_summary:

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


    rule make_qc_plots:

    
        

