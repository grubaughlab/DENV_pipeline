import os
import sys
import datetime as dt
import shutil

from denv_pipeline.utils.misc import *
from denv_pipeline.scripts.make_summary_files import *

denvtype_list = ["DENV1", "DENV2","DENV3","DENV4"]
cwd = os.getcwd()

rule all: #this will be the outputs we want - in the end have the FINAL directory stuff
    input:
        overall_denv_serotype_calls = os.path.join(config["outdir"], "DENV.serotype.calls.tsv")


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
        
        if os.path.exists(os.path.join(params.outdir, "results")):
            if config["overwite"]:
                shutil.rmtree(os.path.join(params.outdir, "results"), ignore_errors=True)
            else:
                sys.stderr.write(f"Error: results file already exists at {os.path.join(params.outdir)}. Use --overwrite flag to delete and regenerate results.")
                sys.exit(-1)

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
                    if sample_dir != "download" and os.path.join(params.indir, sample_dir) != params.tempdir:
                        fw.write(sample_dir + "\n")

rule prepare_jobs:
    output:
        jobs = os.path.join(config["outdir"], "jobs.txt"),
        sample_list = []
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
                    output.sample_list.append(name)
                    basename = name.split("_")[0]
                    primer1 = f"{name}/*R1*"
                    primer2 = f"{name}/*R2*"
                    fw.write(f'bash {input.mapper_script} {basename} {os.path.join(params.indir, primer1)} {os.path.join(params.indir, primer2)} {input.primer_dir} {input.python_script} {params.depth} {params.outdir}')


rule denv_mapper:
    input:
        jobs = rules.prepare_jobs.output.jobs
    output:
        sample_denv_serotype_calls = expand(os.path.join(config["outdir"], "{denv_samples}.serotype.calls.txt"), denv_samples=rules.prepare_jobs.output.sample_list),
        bam_files = expand(os.path.join(config["outdir"], "{denv_type}.{denv_samples}.sort.bam"), denv_type = denvtype_list, denv_samples=rules.prepare_jobs.output.sample_list),
        out_alns = expand(os.path.join(config["outdir"], "{denv_type}.{denv_samples}.20.out.aln"), denv_type = denvtype_list, denv_samples=rules.prepare_jobs.output.sample_list),
        consensus = expand(os.path.join(config["outdir"], "{denv_type}.{denv_samples}.20.cons.fa"), denv_type = denvtype_list, denv_samples=rules.prepare_jobs.output.sample_list)
    run:    
        if config["slurm"]:
            print("preparing for slurm run")
            shell("dsq --job-name denv.mapper --job-file {input.jobs:q} --mem-per-cpu=10G --cpus-per-task=1") 
            filename = f"dsq-jobs-{dt.datetime.today().date()}.sh"
            shell("sbatch {filename}")
        
        else:
            print("running each sample sequentially")
            with open(input.jobs) as f:
                for l in f:
                    command = l.strip("\n")
                    shell("{command}")
            

rule denv_summary:
    input:
        sample_serotype_cals = rules.denv_mapper.output.sample_denv_serotype_calls,
        bam_files = rules.denv_mapper.output.bam_files,
        alignments = rules.denv_mapper.output.out_alns,
        consensus = rules.denv_mapper.output.consensus
    output:
        denv_serotype_calls = os.path.join(config["outdir"], "DENV.serotype.calls.tsv"),
        all_sample_summary = os.path.join(config["outdir"],"summary.all.samples.tsv"),
        top_serotype_calls_all = os.path.join(config["outdir"], "DENV.top.serotype.calls.all.samples.tsv")
    params:
        results_dir = os.path.join(config["outdir"], "results"),
        outdir = config["outdir"],
        python_script = os.path.join(workflow.current_basedir,"make_summary_files.py")
    run:
        os.mkdir(params.results_dir)
        os.mkdir(os.path.join(params.results_dir, "bam_files"))
        os.mkdir(os.path.join(params.results_dir, "variants"))
        os.mkdir(os.path.join(params.results_dir, "depth"))
        os.mkdir(os.path.join(params.results_dir, "consensus_sequences"))

        shell('echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > {output.denv_serotype_calls}')
        shell('cat {params.outdir}/tmp.*.serotype.calls.*.txt >> {output.denv_serotype_calls}')

        make_summary_files.summarise_files(config, os.path.join(params.outdir, "DENV.serotype.calls.tsv"))
        
        shell('echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > {output.all_sample_summary}')
        shell('cat {params.outdir}/*.serotype.calls.txt >> {output.all_sample_summary}')
       
        shell('echo -e "SampleID\tConsSequence\tDepth\tSerotype\tRefSerotypeSequence\tRefSeqLength\tAlignedBases\tCoverageUntrimmed\tCoverageTrimmed" > {output.top_serotype_calls_all}'
        
            'ls {params.outdir}/*.serotype.calls.txt | while read i;'' 
            'do' 
                'cat $i | sort -k8 -n -r | head -1 >> {output.top_serotype_calls_all};' 
            'done'
        )

        


rule tidy_up:
    input:
        serotype_calls = rules.denv_summary.output.denv_serotype_calls,
        all_samples = rules.denv_summary.output.all_sample_summary,
        top_calls_all = rules.denv_summary.output.top_serotype_calls_all
    params:
        temp_files = ["*.cons.qual.txt","*.DENV1.bam", "*.DENV2.bam", "*.DENV3.bam", "*.DENV4.bam", "*.sort.bam.bai", "*.trimmed.bam", "tmp.*.serotype.calls.*.txt",  "*.serotype.txt"],
        tempdir = config["tempdir"],
        resultdir = config["resultdir"]
    run:
        if not config["temp"]:
            for i in temp_files:
                remove_multiple_files(i)
        else:
            for i in temp_files:
                shutil.move(i, {params.tempdir})
        remove_multiple_files("ZZ.tmp000.*")

        shutil.move(input.serotype_calls, {params.resultdir})
        shutil.move(input.all_samples, {params.resultdir})
        shutil.move(input.top_calls_all, {params.resultdir})        


#    # rule make_qc_plots:
#will have to have ct file and column as an argument
# have a secret grubaugh lab option so it picks up all of the samples, not just over 50%

        

