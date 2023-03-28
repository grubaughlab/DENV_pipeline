import os
import sys
import datetime as dt
import shutil

from denv_pipeline.utils.misc import *
from denv_pipeline.scripts.visualisations import *
from denv_pipeline.scripts.make_summary_files import summarise_files


cwd = os.getcwd()

rule all: 
    input:
        os.path.join(config["outdir"], "samples.txt"),
        os.path.join(config["outdir"], "jobs.txt"),
        os.path.join(config["outdir"], "results", "DENV.serotype.calls.tsv"),
        os.path.join(config["outdir"], "results", "variant_plot.pdf")

rule setup:
    output:
        sample_file = os.path.join(config["outdir"], "samples.txt"),
    params:
        tempdir = config["tempdir"],
        indir = config["indir"]
    run:
        with open(output.sample_file, "w") as fw:
            for sample_dir in os.listdir(params.indir):
                if os.path.isdir(os.path.join(params.indir, sample_dir)):
                    if sample_dir != "download" and os.path.join(params.indir, sample_dir) != params.tempdir:
                        fw.write(sample_dir + "\n")

rule prepare_jobs:
    output:
        jobs = os.path.join(config["outdir"], "jobs.txt"),
    input:
        sample_file = rules.setup.output.sample_file,
        mapper_script = os.path.join(workflow.current_basedir,"DENV_MAPPER.sh"),
        primer_dir = config["primer_directory"],
        python_script = os.path.join(workflow.current_basedir,"serotypeCaller.py")
    params:
        outdir = config["outdir"],
        indir = config["indir"],
        depth = config["depth"]
    run:
        with open(output.jobs, 'w') as fw:
            with open(input.sample_file) as f:
                for l in f:
                    name = l.strip("\n")
                    basename = name.split("_")[0]
                    primer1 = f"{name}/*R1*"
                    primer2 = f"{name}/*R2*"
                    fw.write(f'bash {input.mapper_script} {basename} {os.path.join(params.indir, primer1)} {os.path.join(params.indir, primer2)} {input.primer_dir} {input.python_script} {params.depth} {params.outdir}\n')


rule denv_mapper:
    input:
        jobs = rules.prepare_jobs.output.jobs
    output:
        temp_call_files = expand(os.path.join(config["outdir"], "tmp.{sample}.serotype.calls.{depth}.txt"), sample=config["sample_list"], depth=config["depth"]),
        sample_serotype_calls = expand(os.path.join(config["outdir"], "{sample}.serotype.calls.txt"), sample=config["sample_list"]),
        bam_files = expand(os.path.join(config["outdir"], "{sample}.{virus_type}.sort.bam"), sample=config["sample_list"], virus_type=config["option_list"]),
        out_alns = expand(os.path.join(config["outdir"], "{sample}.{virus_type}.{depth}.out.aln"), sample=config["sample_list"], virus_type=config["option_list"], depth=config["depth"]),
        consensus = expand(os.path.join(config["outdir"], "{sample}.{virus_type}.{depth}.cons.fa"), sample=config["sample_list"], virus_type=config["option_list"], depth=config["depth"]),
    params:
        outdir = config["outdir"]
    run:    
        if config["slurm"]:
            print("preparing for slurm run")
            shell("""module load dSQ; 
            dsq --job-name denv.mapper --job-file {input.jobs:q} --mem-per-cpu=10G --cpus-per-task=1""")
             
            filename = f"dsq-jobs-{dt.datetime.today().date()}.sh"
            shell("RES=$(sbatch {filename}) && sbatch --depend=afterok:${RES##* } mapper_done_slurm.sh {params.outdir}") 
        
        else:
            print("running each sample sequentially")
            with open(input.jobs) as f:
                for l in f:
                    command = l.strip("\n")
                    shell("{command}")

rule denv_summary:
    input:
        sample_serotype_calls = rules.denv_mapper.output.sample_serotype_calls,
        bam_files = rules.denv_mapper.output.bam_files,
        alignments = rules.denv_mapper.output.out_alns,
        consensus = rules.denv_mapper.output.consensus,
        temp_call_files = rules.denv_mapper.output.temp_call_files,
    output:
        denv_serotype_calls = os.path.join(config["outdir"], "DENV.serotype.calls.tsv"),
        all_sample_summary = os.path.join(config["outdir"],"summary.all.samples.tsv"),
        top_serotype_calls_all = os.path.join(config["outdir"], "DENV.top.serotype.calls.all.samples.tsv"),
        variant_summary_file = os.path.join(config["outdir"], "results", "variants", "variants_summary.tsv")
    params:
        results_dir = os.path.join(config["outdir"], "results"),
        outdir = config["outdir"],
        python_script = os.path.join(workflow.current_basedir,"make_summary_files.py")
    run:
        make_directory(params.results_dir)
        make_directory(os.path.join(params.results_dir, "bam_files"))
        make_directory(os.path.join(params.results_dir, "variants"))
        make_directory(os.path.join(params.results_dir, "depth"))
        make_directory(os.path.join(params.results_dir, "consensus_sequences"))
        make_directory(os.path.join(params.results_dir, "alignments"))

        shell('echo -e "sample_id\tconsensus_sequence_file\tdepth\tserotype\treference_serotype_name\treference_sequence_length\tnumber_aligned_bases\tcoverage_untrimmed\tcoverage_trimmed" > {output.denv_serotype_calls}')
        shell('cat {input.temp_call_files} >> {output.denv_serotype_calls}')
        
        print("summarising files")
        summarise_files(config, output.denv_serotype_calls)
        
        shell('echo -e "sample_id\tconsensus_sequence_file\tdepth\tserotype\treference_serotype_name\treference_sequence_length\tnumber_aligned_bases\tcoverage_untrimmed\tcoverage_trimmed" > {output.all_sample_summary}')
        shell('cat {input.sample_serotype_calls} >> {output.all_sample_summary}')

        shell('echo -e "sample_id\tconsensus_sequence_file\tdepth\tserotype\treference_serotype_name\treference_sequence_length\tnumber_aligned_bases\tcoverage_untrimmed\tcoverage_trimmed" > {output.top_serotype_calls_all};')
        shell('ls {input.sample_serotype_calls} | while read i; do cat $i | sort -k8 -n -r | head -1 >> {output.top_serotype_calls_all}; done')
        
        alignment_dir = os.path.join(params.results_dir, "alignments")
        for virus_type in config["option_list"]:
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

        serotype_dict, colour_dict, patch_list = prepare_for_plots(input.serotype_calls_file)
        
        variant_plot(params.results_dir, input.variant_summary_file, serotype_dict, colour_dict, patch_list)

        if config["ct_file"] and config["ct_column"] and config["id_column"]:
            ct_plot(params.results_dir, config["ct_file"], config["ct_column"], config["id_column"], input.serotype_calls_file, serotype_dict, colour_dict, patch_list)


rule tidy_up:
    input:
        serotype_calls = rules.denv_summary.output.denv_serotype_calls,
        all_samples = rules.denv_summary.output.all_sample_summary,
        top_calls_all = rules.denv_summary.output.top_serotype_calls_all,
        variant_plot = rules.make_qc_plots.output.variant_plot
    output:
        results_serotype_calls = os.path.join(config["outdir"], "results", "DENV.serotype.calls.tsv")
    params:
        temp_files = ["cons.qual.txt","bam", "sort.bam.bai", "trimmed.bam", "tmp.*.serotype.calls.*.txt", "serotype.calls.txt", "variants.tsv"],
        tempdir = config["tempdir"],
        results_dir = os.path.join(config["outdir"], "results")
    run:
        move_temp_files(config, params.temp_files, params.tempdir)
        clean_up_alignment_components(config, params.tempdir)

        remove_multiple_files("ZZ.tmp000.*")

        shutil.move(input.serotype_calls, params.results_dir)
        shutil.move(input.all_samples, params.results_dir)
        shutil.move(input.top_calls_all, params.results_dir)  

        if config["download"]:
            make_directory(os.path.join(config["outdir"], "downloads")) 
            for directory in os.listdir(params.results_dir):
                if directory != "bam_files":
                    source = os.path.join(params.results_dir, directory)
                    dest = os.path.join(config["outdir"], "downloads")
                    shell("cp -r {source} {dest}")


        

