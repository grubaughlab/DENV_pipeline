import os
import sys
import datetime as dt
import shutil

from denv_pipeline.utils import misc
from denv_pipeline.scripts import visualisations
from denv_pipeline.scripts import make_summary_files


cwd = os.getcwd()

rule all: 
    input:
        os.path.join(config["outdir"], "results", "virus_calls.tsv"),
        os.path.join(config["outdir"], "results", "variant_plot.pdf")

rule mapper:
    input:
        read_location = os.path.join(config["indir"], "{sample}")
    output:
        individual_all_virustype_info = temp(os.path.join(config["tempdir"], "{sample}_all_virustype_info.txt"))
    log:
        log = os.path.join(config["outdir"], "log_files", "_".join(["{sample}", "mapping.log"]))
    params:
        mapper_script = os.path.join(workflow.current_basedir,"mapper.sh"),
        primer_dir = config["reference_directory"],
        depth = config["depth"],
        tempdir = config["tempdir"],
        python_script = os.path.join(workflow.current_basedir,"serotype_caller.py"),
        python_script2 = os.path.join(workflow.current_basedir, "make_empty_files.py")
    resources:
        partition="day",
        mem_mb_per_cpu="10G",
        cpus_per_task=2,
        runtime=300
    run:
        shell("{params.mapper_script} {wildcards.sample} {input.read_location}/*R1* {input.read_location}/*R2* {params.primer_dir} {params.python_script} {params.python_script2} {params.depth} {params.tempdir} {log.log}  >> {log.log} 2>&1")
        
        if not os.path.exists(os.path.join(params.tempdir,f"{wildcards.sample}_all_virustype_info.txt")):
            shell("touch {params.tempdir}/{wildcards.sample}_all_virustype_info.txt")

rule summary:
    input:
    #have to be like this (ie not rules.output) otherwise the wildcards don't work
        individual_all_virustype_info = expand(os.path.join(config["tempdir"], "{sample}_all_virustype_info.txt"), sample=config["sample_list"])
    output:
        serotype_calls = os.path.join(config["outdir"], "results", "virus_calls.tsv"),
        all_serotype_summary = os.path.join(config["outdir"], "results", "summary_all_samples.tsv"),
        top_serotype_summary = os.path.join(config["outdir"], "results", "top_virus_all_samples.tsv"),
        variant_summary_file = os.path.join(config["outdir"], "results", "variants", "variants_summary.tsv")
    params:
        results_dir = os.path.join(config["outdir"], "results"),
        outdir = config["outdir"],
        tempdir = config["tempdir"],
        python_script = os.path.join(workflow.current_basedir,"make_summary_files.py")
    run:        
        print("summarising files")
        make_summary_files.summarise_files(config, input.individual_all_virustype_info, output.serotype_calls, output.top_serotype_summary, output.all_serotype_summary)
        
        if config["download"]:
            for directory in os.listdir(params.results_dir):
                if directory != "bam_files":
                    source = os.path.join(params.results_dir, directory)
                    dest = os.path.join(config["outdir"], "downloads")
                    shell("cp -r {source} {dest}")

        if not config["temp"]:
            shutil.rmtree(config["tempdir"])

rule make_qc_plots:
    input:
        serotype_calls_file = rules.summary.output.serotype_calls,
        variant_summary_file = rules.summary.output.variant_summary_file
    output:
        variant_plot = os.path.join(config["outdir"], "results", "variant_plot.pdf")
    params:
        results_dir = rules.summary.params.results_dir
    run:
        serotype_dict, colour_dict, patch_list = visualisations.prepare_for_plots(input.serotype_calls_file)
        
        visualisations.variant_plot(params.results_dir, input.variant_summary_file, serotype_dict, colour_dict, patch_list)

        if config["ct_file"] and config["ct_column"] and config["id_column"]:
            visualisations.ct_plot(params.results_dir, config["ct_file"], config["ct_column"], config["id_column"], input.serotype_calls_file, serotype_dict, colour_dict, patch_list)


        

