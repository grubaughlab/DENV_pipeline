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
        os.path.join(config["outdir"], "results", "DENV.serotype.calls.tsv"),
        os.path.join(config["outdir"], "results", "variant_plot.pdf")

rule denv_mapper:
    input:
        read_location = os.path.join(config["indir"], "{sample}")
    output:
        individual_all_virustype_info = temp(os.path.join(config["outdir"], "{sample}_all_virustype_info.txt"))
    log:
        log = os.path.join(config["outdir"], "log_files", "_".join(["{sample}", "mapping.log"]))
    params:
        mapper_script = os.path.join(workflow.current_basedir,"DENV_MAPPER.sh"),
        primer_dir = config["reference_directory"],
        depth = config["depth"],
        outdir = config["outdir"],
        python_script = os.path.join(workflow.current_basedir,"serotypeCaller.py")
    resources:
        partition="general",
        mem_mb_per_cpu="10G",
        cpus_per_task=1,
        runtime=300
    run:
        shell("{params.mapper_script} {wildcards.sample} {input.read_location}/*R1* {input.read_location}/*R2* {params.primer_dir} {params.python_script} {params.depth} {params.outdir} {log.log}  >> {log.log} 2>&1")
        
        if not os.path.exists(os.path.join(params.outdir,f"{wildcards.sample}_all_virustype_info.txt")):
            shell("touch {params.outdir}/{wildcards.sample}_all_virustype_info.txt")

rule denv_summary:
    input:
    #have to be like this (ie not rules.output) otherwise the wildcards don't work
        individual_all_virustype_info = expand(os.path.join(config["outdir"], "{sample}_all_virustype_info.txt"), sample=config["sample_list"])
    output:
        serotype_calls = os.path.join(config["outdir"], "results", "DENV.serotype.calls.tsv"),
        all_serotype_summary = os.path.join(config["outdir"], "results", "summary.all.samples.tsv"),
        top_serotype_summary = os.path.join(config["outdir"], "results", "DENV.top.serotype.calls.all.samples.tsv"),
        variant_summary_file = os.path.join(config["outdir"], "results", "variants", "variants_summary.tsv")
    params:
        results_dir = os.path.join(config["outdir"], "results"),
        outdir = config["outdir"],
        python_script = os.path.join(workflow.current_basedir,"make_summary_files.py")
    log:
        os.path.join(config["outdir"], "log_files", "summary.log")
    run:        
        print("summarising files")
        make_summary_files.summarise_files(config, input.individual_all_virustype_info, output.serotype_calls, output.top_serotype_summary, output.all_serotype_summary)
        
        alignment_dir = os.path.join(params.results_dir, "alignments")
        for virus_type in config["virus_type_list"]:
            for file_name in os.listdir(alignment_dir):
                if virus_type in file_name:
                    if "trim" in file_name:
                        shell("cat {alignment_dir}/{file_name} >> {alignment_dir}/{virus_type}.trim.aln")
                    else:
                        shell("cat {alignment_dir}/{file_name} >> {alignment_dir}/{virus_type}.untrim.aln")

rule make_qc_plots:
    input:
        serotype_calls_file = rules.denv_summary.output.serotype_calls,
        variant_summary_file = rules.denv_summary.output.variant_summary_file
    output:
        variant_plot = os.path.join(config["outdir"], "results", "variant_plot.pdf")
    params:
        results_dir = rules.denv_summary.params.results_dir
    run:
        serotype_dict, colour_dict, patch_list = visualisations.prepare_for_plots(input.serotype_calls_file)
        
        visualisations.variant_plot(params.results_dir, input.variant_summary_file, serotype_dict, colour_dict, patch_list)

        if config["ct_file"] and config["ct_column"] and config["id_column"]:
            visualisations.ct_plot(params.results_dir, config["ct_file"], config["ct_column"], config["id_column"], input.serotype_calls_file, serotype_dict, colour_dict, patch_list)


rule tidy_up:
    input:
        serotype_calls = rules.denv_summary.output.serotype_calls,
        all_serotype_summary = rules.denv_summary.output.all_serotype_summary,
        top_serotype_summary = rules.denv_summary.output.top_serotype_summary,
        variant_plot = rules.make_qc_plots.output.variant_plot
    output:
        results_serotype_calls = os.path.join(config["outdir"], "results", "DENV.serotype.calls.tsv")
    params:
        temp_files = ["cons.qual.txt", "bam", "bam.bai", "trimmed.bam", "tmp.*.serotype.calls.*.txt", "serotype.calls.txt", "variants.tsv"],
        tempdir = config["tempdir"],
        results_dir = os.path.join(config["outdir"], "results")
    run:
        make_summary_files.move_temp_files(config, params.temp_files, params.tempdir)
        make_summary_files.clean_up_alignment_components(config, params.tempdir)
        misc.remove_multiple_files(os.path.join(config["outdir"], "ZZ.tmp000.*"))

        if config["download"]:
            for directory in os.listdir(params.results_dir):
                if directory != "bam_files":
                    source = os.path.join(params.results_dir, directory)
                    dest = os.path.join(config["outdir"], "downloads")
                    shell("cp -r {source} {dest}")


        

