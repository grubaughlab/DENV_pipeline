import os
import sys

rule all: #this will be the outputs we want
    input:
        os.path.join(config["cwd"], "jobs.txt")

rule setup:
    output:
        sample_file = os.path.join(config["cwd"], "samples.txt")
    params:
        cwd = config["cwd"]
    run:
        if not os.path.exists(config["cwd"]):
            os.mkdir(config["cwd"])
        
        shell("cd {params.cwd}")
        if os.path.exists("DENV.serotype.calls.tsv"):
            os.remove("DENV.serotype.calls.tsv")
        
        
        shell("/home/bioinfo/software/knightlab/bin_Mar2018/ycgaFastq {config['symlink']}")
        shell("ls | grep -v samples > {output.sample_file:q}")


    
rule denv_mapper:
    input:
        mapper_script = os.path.join(workflow.current_basedir,"DENV_MAPPER.sh"),
        sample_file = rules.setup.output.sample_file,
        refs = os.path.join(config["denv_primers"], "DENV.refs.txt")
    output:
        jobs = os.path.join(config["cwd"], "jobs.txt")
        #have the full denv ones here
    run:
        with open(output.jobs, 'w') as fw:
            with open(input.sample_file) as f:
                for l in f:
                    name = l.strip("\n")
                    fw.write(f"bash {input.mapper_script} {name}/*/{name}*_R1_*.fastq.gz {name}/*/{name}*_R2_*.fastq.gz {input.refs}")

            
        if config["slurm"]:
            print("preparing for slurm run")
            shell("dsq --job-name denv.mapper --job-file output.jobs:q --mem-per-cpu=10G --cpus-per-task=1") 
            filename = shell("ls | grep dsq")
            shell("sbatch filename")
        
        else:
            print("running each sample sequentially")
            with open(output.jobs) as f:
                for l in f:
                    command = l.strip("\n")
                    shell("{command}")
        
            



