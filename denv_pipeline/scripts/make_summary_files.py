import os
import csv
import argparse
from collections import defaultdict
import shutil

def summarise_files(config, serotype_calls):

    min_coverage = []
    serotypes = defaultdict(list)
    with open(serotype_calls) as f:
        data = csv.DictReader(f, delimiter="\t")
        for l in data:
            if float(l['coverage_untrimmed']) >= 50:
                min_coverage.append(l)
                serotypes[l['sample_id']].append(l['serotype'])
    
    sort_variant_files(config, serotypes)
    get_right_serotype_files(config, serotypes)
    
    if config["temp"]:
        with open(os.path.join(config["tempdir"], "DENV.serotype.calls.mincov50.final.tsv"), 'w') as fw:
            headers = ["sample_id", "consensus_sequence_file", "depth", "serotype", "reference_serotype_name", "reference_sequence_length", "number_aligned_bases", "coverage_untrimmed", "coverage_trimmed"]
            writer = csv.DictWriter(fw, headers, delimiter="\t")
            writer.writeheader()
            write_dict = {}
            for i in min_coverage:
                for k,v in i.items():
                    write_dict[k] = v
                writer.writerow(write_dict)

    return


def sort_variant_files(config, serotypes):
    
    old_to_new = {'POS': 'position', 'REF': 'reference_base', 'ALT': 'alternative_base', 'REF_DP': 'reference_depth', 'REF_RV': 'reference_depth_reverse', 'REF_QUAL': 'reference_quality', 'ALT_DP': 'alternate_depth', 'ALT_RV': 'alternate_depth_reverse', 'ALT_QUAL': 'alternative_quality', 'ALT_FREQ': 'alternative_frequency', 'TOTAL_DP': 'total_depth', 'PVAL': 'p_value_fisher', 'PASS': 'pass', 'GFF_FEATURE': 'gff_feature', 'REF_CODON': 'reference_codon', 'REF_AA': 'reference_amino_acid', 'ALT_CODON': 'alternative_codon', 'ALT_AA': 'alternative_amino_acid'}

    summary_file = open(os.path.join(config["outdir"], "results", "variants", "variants_summary.tsv"), 'w')
    summary_file.write("sample_id\tserotype\tvariant_count\n")

    for file in os.listdir(config["outdir"]):
        if file.endswith(".variants.tsv"):
            count = 0
            new_file_name = file.replace(".variants.tsv", ".variants_frequency.tsv")
            new_file = os.path.join(config['outdir'], new_file_name)
            with open(new_file, 'w') as fw:
                headers = list(old_to_new.values())
                writer = csv.DictWriter(fw, headers, delimiter="\t")
                writer.writeheader()
                with open(os.path.join(config['outdir'], file)) as f:
                    data = csv.DictReader(f, delimiter="\t")
                    for l in data: 
                        write_dict = {}
                        if l['PASS'] == "TRUE" and float(l['ALT_FREQ']) < 0.8 and float(l['ALT_FREQ']) > 0.2:
                            for old,new in old_to_new.items():
                                write_dict[new] = l[old]
                            writer.writerow(write_dict)
                            count += 1

            sero = file.split(".")[1]
            sample = file.split(".")[0]
            if sero in serotypes[sample]:
                summary_file.write(f"{sample}\t{sero}\t{count}\n")

    summary_file.close()
    
    return

def clean_depth_file(config, depth_file):

    new_depth_file = depth_file.split("/")[-1]
    headers = ["position", "depth"]
    with open(new_depth_file, 'w') as fw:
        writer = csv.DictWriter(fw,headers, delimiter="\t")
        writer.writeheader()
        with open(depth_file) as f:
            for l in f:
                write_dict = {}
                toks = l.strip("\n").split("\t")
                write_dict["position"] = toks[1]
                write_dict["depth"] = toks[2]
                writer.writerow(write_dict)

    return new_depth_file

    
def get_right_serotype_files(config, serotypes):

    bam_files = set()
    consensus = set()
    depths = set()
    variant_frequencies = set()    
    alignments = set()
    depth = config["depth"]
    full_virus_type_list = config["virus_type_list"]
    
    unwanted = []
    for sample, serotype_lst in serotypes.items():
        for option in full_virus_type_list:
            
            bam_file = f'{sample}.{option}.sort.bam'
            consensus_file = f'{sample}.{option}.{depth}.cons.fa'
            depth_file = f'{sample}.{option}.depth.txt'
            variant_frequency = f'{sample}.{option}.{depth}.variants_frequency.tsv'
            trimmed = f'{sample}.{option}.{depth}.out.trim.aln'
            untrimmed = f'{sample}.{option}.{depth}.out.aln'

            if option in serotype_lst:
                bam_files.add(bam_file)
                consensus.add(consensus_file)
                depths.add(depth_file)
                variant_frequencies.add(variant_frequency)
                alignments.add(trimmed)
                alignments.add(untrimmed)
            else:
                unwanted.append(bam_file)
                unwanted.append(consensus_file)
                unwanted.append(depth_file)
                unwanted.append(variant_frequency)
                unwanted.append(trimmed)
                unwanted.append(untrimmed)

    for bam in bam_files:
        source = os.path.join(config['outdir'], bam)
        dest = os.path.join(config["outdir"], "results", "bam_files")
        shutil.move(source, dest)

    for cons in consensus:
        source = os.path.join(config['outdir'], cons)
        dest = os.path.join(config["outdir"], "results", "consensus_sequences")
        shutil.move(source, dest)

    for dep in depths:
        source = clean_depth_file(config, os.path.join(config['outdir'], dep))
        dest = os.path.join(config["outdir"], "results", "depth")
        shutil.move(source, dest)
        unwanted.append(dep)

    for var_freq in variant_frequencies:
        source = os.path.join(config['outdir'], var_freq)
        dest = os.path.join(config["outdir"], "results", "variants")
        shutil.move(source, dest)

    for aln in alignments:
        source = os.path.join(config['outdir'], aln)
        dest = os.path.join(config["outdir"], "results", "alignments")
        shutil.move(source, dest)

    for i in unwanted:
        if config["temp"]:
            shutil.move(os.path.join(config['outdir'], i), config['tempdir'])
        else:
            os.remove(os.path.join(config['outdir'], i))


def move_temp_files(config, temp_files, dest):

    contains_depth = ["cons.qual.txt", "variants.tsv"]
    depth = config["depth"]
    for file_pattern in temp_files:
        for sample in config["sample_list"]:
            if "serotype.calls" in file_pattern:
                if "tmp" in file_pattern:
                    name = f"tmp.{sample}.serotype.calls.{depth}.txt"
                else:
                    name = f"{sample}.serotype.calls.txt"

                if config["temp"]:
                    shutil.move(os.path.join(config["outdir"], name), dest)
                else:
                    os.remove(os.path.join(config["outdir"], name))
            
            else:
                for option in config["virus_type_list"]:
                    if ".bam" in file_pattern:
                        if ".sort" in file_pattern and ".bai" not in file_pattern:
                            pass
                        else:
                            name = f"{sample}.{option}.{file_pattern}"
                    else:
                        if file_pattern in contains_depth:
                            name = f"{sample}.{option}.{depth}.{file_pattern}"
                        else:
                            name = f"{sample}.{option}.{file_pattern}"

                    if config["temp"]:
                        if "bam.bai" in file_pattern:
                            if os.path.exists(os.path.join(config["outdir"], name)):
                                shutil.move(os.path.join(config["outdir"], name), dest)
                            else:
                                pass
                        else:
                            shutil.move(os.path.join(config["outdir"], name), dest)
                    else:
                        if "bam.bai" in file_pattern:
                            if os.path.exists(os.path.join(config["outdir"], name)):
                                os.remove(os.path.join(config["outdir"], name))
                            else:
                                pass
                        else:
                            os.remove(os.path.join(config["outdir"], name))

    

def clean_up_alignment_components(config, temp_dir):

    path = os.path.join(config["outdir"],"results", "alignments")
    depth = config["depth"]
    components = []
    for option in config["virus_type_list"]:
        for sample in config["sample_list"]:
            file1 = f'{sample}.{option}.{depth}.out.aln'
            file2 = f'{sample}.{option}.{depth}.out.trim.aln'

            if os.path.exists(os.path.join(path, file1)):
                if config["temp"]:
                    shutil.move(os.path.join(path, file1), temp_dir)
                    shutil.move(os.path.join(path, file2), temp_dir)
                else:
                    os.remove(os.path.join(path, file1))
                    os.remove(os.path.join(path, file2))

    #if there's no trimmed bedfile then delete the empty trimmed alignment
    for alignment in path:
        if alignment.endswith(".trim.aln") and os.stat(os.path.join(path, alignment)).st_size == 0:
            if config["temp"]:
                shutil.move(os.path.join(path, alignment), temp_dir)
            else:
                os.remove(os.path.join(path, alignment))