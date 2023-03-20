import os
import sys
import datetime as dt

from denv_pipeline.utils.misc import *

denvtype_list = ["DENV1", "DENV2","DENV3","DENV4"]
cwd = os.getcwd()

# rule all: #this will be the outputs we want - in the end have the FINAL directory stuff
#     input:
#         overall_denv_serotype_calls = os.path.join(config["outdir"], "DENV.serotype.calls.tsv")


rule setup:
    output:
        sample_file = os.path.join(config["outdir"], "samples.txt"),
    params:
        outdir = config["outdir"],
        symlink = config["symlink"],
        tempdir = config["tempdir"],
        indir = config["indir"]
    run:
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        if config["download"] and not os.path.exists(os.path.join(params.outdir, "downloads")):
           os.mkdir(os.path.join(params.outdir, "downloads"))
        if config["temp"] and not os.path.exists(params.tempdir):
           os.mkdir(params.tempdir)

        if params.symlink != "":
            shell("cd {params.indir}")
            shell("/home/bioinfo/software/knightlab/bin_Mar2018/ycgaFastq {params.symlink:q}")
            for sample_dir in os.listdir(params.indir):
                if os.path.isdir(sample_dir):
                    folder_path = os.path.join(sample_dir, "Unaligned")
                    shell("mv {folder_path}/* {sample_dir}/")
                    shell("rm -r {folder_path}")
            shell("cd {cwd}")

        with open(output.sample_file, "w") as fw:
            for sample_dir in os.listdir(params.indir):
                if os.path.isdir(os.path.join(params.indir, sample_dir)):
                    fw.write(sample_dir + "\n")

# rule prepare_jobs:
#     output:
#         jobs = os.path.join(config["outdir"], "jobs.txt"),
#         denv_sample_list = []
#     input:
#         sample_file = rules.setup.output.sample_file,
#         mapper_script = os.path.join(workflow.current_basedir,"DENV_MAPPER.sh"),
#         primer_dir = config["denv_primers"],
#         python_script = os.path.join(workflow.current_basedir,"serotypeCaller.py"),
#         outdir = config["outdir"],
#         indir = config["indir"]
#     run:
#         with open(output.jobs, 'w') as fw:
#             with open(input.sample_file) as f:
#                 for l in f:
#                     name = l.strip("\n")
#                     output.denv_sample_list.append(name)
#                     basename = name.split("_")[0]
#                     primer1 = f"{name}/*R1*"
#                     primer2 = f"{name}/*R2*"
#                     fw.write(f'bash {input.mapper_script} {basename} {os.path.join(input.indir, primer1)} {os.path.join(input.indir, primer2)} {input.primer_dir} {input.python_script} {input.outdir}')


# rule denv_mapper:
#     input:
#         jobs = rules.prepare_jobs.output.jobs
#     output:
#         overall_denv_serotype_calls = os.path.join(config["outdir"], "DENV.serotype.calls.tsv"),
#         sample_denv_serotype_calls = expand(os.path.join(config["outdir"], "{denv_samples}.serotype.calls.txt"), denv_samples=rules.prepare_jobs.output.denv_sample_list),
#         frequency_counts = expand(os.path.join(config["outdir"], "{denv_type}.{denv_samples}.20_variants_frequency_count.txt"), denv_type = denvtype_list, denv_samples=rules.prepare_jobs.output.denv_sample_list),
#         bam_files = expand(os.path.join(config["outdir"], "{denv_type}.{denv_samples}.sort.bam"), denv_type = denvtype_list, denv_samples=rules.prepare_jobs.output.denv_sample_list),
#         out_alns = expand(os.path.join(config["outdir"], "{denv_type}.{denv_samples}.20.out.aln"), denv_type = denvtype_list, denv_samples=rules.prepare_jobs.output.denv_sample_list),
#         consensus = expand(os.path.join(config["outdir"], "{denv_type}.{denv_samples}.20.cons.fa"), denv_type = denvtype_list, denv_samples=rules.prepare_jobs.output.denv_sample_list)
        
#     run:    
#         if config["slurm"]:
#             print("preparing for slurm run")
#             shell("dsq --job-name denv.mapper --job-file {input.jobs:q} --mem-per-cpu=10G --cpus-per-task=1") 
#             filename = f"dsq-jobs-{dt.datetime.today().date()}.sh"
#             shell("sbatch {filename}")
        
#         else:
#             print("running each sample sequentially")
#             with open(input.jobs) as f:
#                 for l in f:
#                     command = l.strip("\n")
#                     shell("{command}")
            

#    # rule denv_summary:
# #
#  #       input:

#   #      output:

#    #     parameters:


#     #    run:
#         #DENV_summarise.sh


    
# rule prepare_outputs:
#     params:
#         temp_files = ["*.cons.qual.txt","*.DENV1.bam", "*.DENV2.bam", "*.DENV3.bam", "*.DENV4.bam", "*.sort.bam.bai", "*.trimmed.bam", "tmp.*.serotype.calls.*.txt",  "*.serotype.txt"],
#         tempdir = config["tempdir"]
#     run:
#         if not config["temp"]:
#             for i in temp_files:
#                 remove_multiple_files(i)
#         else:
#             for i in temp_files:
#                 shell("mv {i} {params.tempdir}")
#         remove_multiple_files("ZZ.tmp000.*")


#         #make FINAL directory, call it "results" 
        


#    # rule make_qc_plots:


        

