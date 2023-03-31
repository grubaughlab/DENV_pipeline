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
        primer1 = (os.path.join(config["indir"], "{sample}", "".join(["{sample}", config["fastq_filestem_R1"]]))),
        primer2 = (os.path.join(config["indir"], "{sample}", "".join(["{sample}", config["fastq_filestem_R1"]])))
    output:
        temp_call_files = (os.path.join(config["outdir"], ".".join(["tmp.{sample}.serotype.calls", str(config["depth"]), "txt"]))),
        sample_serotype_calls = (os.path.join(config["outdir"], "{sample}.serotype.calls.txt"))
    log:
        log = os.path.join(config["outdir"], "log_files", "_".join(["{sample}", "mapping.log"]))
    params:
        mapper_script = os.path.join(workflow.current_basedir,"DENV_MAPPER.sh"),
        primer_dir = config["primer_directory"],
        depth = config["depth"],
        outdir = config["outdir"],
        python_script = os.path.join(workflow.current_basedir,"serotypeCaller.py"),
    resources:
        partition="general",
        mem_mb_per_cpu="10G",
        cpus_per_task=1
    shell:
        "{params.mapper_script} {wildcards.sample} {input.primer1} {input.primer2} {params.primer_dir} {params.python_script} {params.depth} {params.outdir} {log.log}  >> {log.log} 2>&1"

rule denv_summary:
    input:
        sample_serotype_calls = expand(os.path.join(config["outdir"], "{sample}.serotype.calls.txt"), sample=config["sample_list"]),
        temp_call_files = expand(os.path.join(config["outdir"], "tmp.{sample}.serotype.calls.{depth}.txt"), sample=config["sample_list"], depth=config["depth"])
    output:
        denv_serotype_calls = os.path.join(config["outdir"], "DENV.serotype.calls.tsv"),
        all_sample_summary = os.path.join(config["outdir"],"summary.all.samples.tsv"),
        top_serotype_calls_all = os.path.join(config["outdir"], "DENV.top.serotype.calls.all.samples.tsv"),
        variant_summary_file = os.path.join(config["outdir"], "results", "variants", "variants_summary.tsv")
    params:
        results_dir = os.path.join(config["outdir"], "results"),
        outdir = config["outdir"],
        python_script = os.path.join(workflow.current_basedir,"make_summary_files.py")
    log:
        os.path.join(config["outdir"], "log_files", "summary.log")
    run:
        misc.make_directory(params.results_dir)
        misc.make_directory(os.path.join(params.results_dir, "bam_files"))
        misc.make_directory(os.path.join(params.results_dir, "variants"))
        misc.make_directory(os.path.join(params.results_dir, "depth"))
        misc.make_directory(os.path.join(params.results_dir, "consensus_sequences"))
        misc.make_directory(os.path.join(params.results_dir, "alignments"))

        shell('echo -e "sample_id\tconsensus_sequence_file\tdepth\tserotype\treference_serotype_name\treference_sequence_length\tnumber_aligned_bases\tcoverage_untrimmed\tcoverage_trimmed" > {output.denv_serotype_calls}')
        shell('cat {input.temp_call_files} >> {output.denv_serotype_calls}')
        
        print("summarising files")
        make_summary_files.summarise_files(config, output.denv_serotype_calls)
        
        shell('echo -e "sample_id\tconsensus_sequence_file\tdepth\tserotype\treference_serotype_name\treference_sequence_length\tnumber_aligned_bases\tcoverage_untrimmed\tcoverage_trimmed" > {output.all_sample_summary}')
        shell('cat {input.sample_serotype_calls} >> {output.all_sample_summary}')

        shell('echo -e "sample_id\tconsensus_sequence_file\tdepth\tserotype\treference_serotype_name\treference_sequence_length\tnumber_aligned_bases\tcoverage_untrimmed\tcoverage_trimmed" > {output.top_serotype_calls_all};')
        shell('ls {input.sample_serotype_calls} | while read i; do cat $i | sort -k8 -n -r | head -1 >> {output.top_serotype_calls_all}; done')
        
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
        serotype_calls_file = rules.denv_summary.output.denv_serotype_calls,
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
        serotype_calls = rules.denv_summary.output.denv_serotype_calls,
        all_samples = rules.denv_summary.output.all_sample_summary,
        top_calls_all = rules.denv_summary.output.top_serotype_calls_all,
        variant_plot = rules.make_qc_plots.output.variant_plot
    output:
        results_serotype_calls = os.path.join(config["outdir"], "results", "DENV.serotype.calls.tsv")
    params:
        temp_files = ["cons.qual.txt","bam", "sort.bam.bai", "bam.bai", "trimmed.bam", "tmp.*.serotype.calls.*.txt", "serotype.calls.txt", "variants.tsv"],
        tempdir = config["tempdir"],
        results_dir = os.path.join(config["outdir"], "results")
    run:
        make_summary_files.move_temp_files(config, params.temp_files, params.tempdir)
        make_summary_files.clean_up_alignment_components(config, params.tempdir)

        misc.remove_multiple_files(os.path.join(config["outdir"], "ZZ.tmp000.*"))

        shutil.move(input.serotype_calls, params.results_dir)
        shutil.move(input.all_samples, params.results_dir)
        shutil.move(input.top_calls_all, params.results_dir)  

        if config["download"]:
            misc.make_directory(os.path.join(config["outdir"], "downloads")) 
            for directory in os.listdir(params.results_dir):
                if directory != "bam_files":
                    source = os.path.join(params.results_dir, directory)
                    dest = os.path.join(config["outdir"], "downloads")
                    shell("cp -r {source} {dest}")


        

