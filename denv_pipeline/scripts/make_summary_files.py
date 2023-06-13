import os
import csv
import argparse
from collections import defaultdict
import shutil

def summarise_files(config, per_sample_files, serotype_call_file, top_call_file, all_info_file):

    serotypes = defaultdict(list)
    all_lines = []
    top_calls = []
    serotype_call = []

    for file in per_sample_files:
        found_top = False
        possible_tops = []
        with open(file) as f:
            data = csv.DictReader(f, delimiter="\t")
            for l in data:
                possible_tops.append(l)
                all_lines.append(l)

                if l['coverage_untrimmed'] != "0":
                    if float(l['coverage_untrimmed']) >= 50:
                        serotype_call.append(l)
                        top_calls.append(l)
                        found_top = True
                        serotypes[l['sample_id']] = l['serotype']
            
        if not found_top:
            if len(possible_tops) > 0:
                top = sorted(possible_tops, key=lambda x:float(x['coverage_untrimmed']), reverse=True)[0]
                top_calls.append(top)
            else:
                continue

    headers = ["sample_id","consensus_sequence_file","depth","serotype","reference_sequence_name","reference_sequence_length","number_aligned_bases","coverage_untrimmed","coverage_trimmed"]
    
    with open(serotype_call_file, 'w') as fw:
        writer = csv.DictWriter(fw,delimiter="\t", fieldnames=headers)
        writer.writeheader()
        for line in serotype_call:
            writer.writerow(line)

    with open(top_call_file, 'w') as fw:
        writer = csv.DictWriter(fw,delimiter="\t", fieldnames=headers)
        writer.writeheader()
        for line in top_calls:
            writer.writerow(line)

    with open(all_info_file, 'w') as fw:
        writer = csv.DictWriter(fw,delimiter="\t", fieldnames=headers)
        writer.writeheader()
        for line in all_lines:
            writer.writerow(line)


    sort_variant_files(config, serotypes)
    alignments = get_right_serotype_files(config, serotypes)
    make_alignments(config, alignments)

    return


def sort_variant_files(config, serotypes):

    old_to_new = {'POS': 'position', 'REF': 'reference_base', 'ALT': 'alternative_base', 'REF_DP': 'reference_depth', 'REF_RV': 'reference_depth_reverse', 'REF_QUAL': 'reference_quality', 'ALT_DP': 'alternate_depth', 'ALT_RV': 'alternate_depth_reverse', 'ALT_QUAL': 'alternative_quality', 'ALT_FREQ': 'alternative_frequency', 'TOTAL_DP': 'total_depth', 'PVAL': 'p_value_fisher', 'PASS': 'pass', 'GFF_FEATURE': 'gff_feature', 'REF_CODON': 'reference_codon', 'REF_AA': 'reference_amino_acid', 'ALT_CODON': 'alternative_codon', 'ALT_AA': 'alternative_amino_acid'}

    summary_file = open(os.path.join(config["outdir"], "results", "variants", "variants_summary.tsv"), 'w')
    summary_file.write("sample_id\tserotype\tvariant_count\n")

    for file in os.listdir(config["tempdir"]):
        if file.endswith(".variants.tsv"):
            count = 0
            new_file_name = file.replace(".variants.tsv", ".variants_frequency.tsv")
            new_file = os.path.join(config['tempdir'], new_file_name)
            with open(new_file, 'w') as fw:
                headers = list(old_to_new.values())
                writer = csv.DictWriter(fw, headers, delimiter="\t")
                writer.writeheader()
                with open(os.path.join(config['tempdir'], file)) as f:
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
    bam_indices = set()
    consensus = set()
    depths = set()
    variant_frequencies = set()    
    alignments = set()
    depth = config["depth"]
    full_virus_type_list = config["virus_type_list"]
    
    for sample, serotype_lst in serotypes.items():
        for option in full_virus_type_list:
            
            bam_file = f'{sample}.{option}.sort.bam'
            bam_index = f'{sample}.{option}.sort.bam.bai'
            consensus_file = f'{sample}.{option}.{depth}.cons.fa'
            depth_file = f'{sample}.{option}.depth.txt'
            variant_frequency = f'{sample}.{option}.{depth}.variants_frequency.tsv'
            trimmed = f'{sample}.{option}.{depth}.out.trim.aln'
            untrimmed = f'{sample}.{option}.{depth}.out.aln'

            if option in serotype_lst:
                bam_files.add(bam_file)
                bam_indices.add(bam_index)
                consensus.add(consensus_file)
                depths.add(depth_file)
                variant_frequencies.add(variant_frequency)
                alignments.add(trimmed)
                alignments.add(untrimmed)

    for bam in bam_files:
        source = os.path.join(config['tempdir'], bam)
        dest = os.path.join(config["outdir"], "results", "bam_files")
        shutil.move(source, dest)
    
    for bam in bam_indices:
        source = os.path.join(config['tempdir'], bam)
        dest = os.path.join(config["outdir"], "results", "bam_files")
        shutil.move(source, dest)

    for cons in consensus:
        source = os.path.join(config['tempdir'], cons)
        dest = os.path.join(config["outdir"], "results", "consensus_sequences")
        shutil.move(source, dest)

    for dep in depths:
        source = clean_depth_file(config, os.path.join(config['tempdir'], dep))
        dest = os.path.join(config["outdir"], "results", "depth")
        shutil.move(source, dest)

    for var_freq in variant_frequencies:
        source = os.path.join(config['tempdir'], var_freq)
        dest = os.path.join(config["outdir"], "results", "variants")
        shutil.move(source, dest)

    return alignments


def make_alignments(config, alignments):

    alignment_dir = os.path.join(config["outdir"], "results", "alignments")
    tempdir = config["tempdir"]

    for virus_type in config["virus_type_list"]:
        new_trimmed_file = f'{alignment_dir}/{virus_type}.trim.aln'
        new_untrimmed_file = f'{alignment_dir}/{virus_type}.untrim.aln'
        for aln in alignments:
            if not os.stat(os.path.join(tempdir, aln)).st_size == 0:
                if virus_type in aln:
                    if "trim" in aln:
                        os.system(f"cat {tempdir}/{aln} >> {new_trimmed_file}")
                    else:
                        os.system(f"cat {tempdir}/{aln} >> {new_untrimmed_file}")


def make_empty_file(tempdir, depth, sample_name, serotype):
    
    headers = ["sample_id","consensus_sequence_file","depth","serotype","reference_sequence_name","reference_sequence_length","number_aligned_bases","coverage_untrimmed","coverage_trimmed"]

    outfile = os.path.join(tempdir, f'{sample_name}_all_virustype_info.txt')

    write_dict = {}
    write_dict["sample_id"] = sample_name
    write_dict["consensus_sequence_file"] = "NA"
    write_dict['depth'] = depth

    write_dict["serotype"] = "NA"
    write_dict["reference_sequence_name"] = serotype
    write_dict["reference_sequence_length"] = "NA"
    write_dict["number_aligned_bases"] = 0
    write_dict["coverage_untrimmed"] = 0
    write_dict["coverage_trimmed"] = 0

    if not os.path.exists(outfile):
        with open(outfile, 'w') as fw:
            writer = csv.DictWriter(fw, delimiter = "\t", fieldnames=headers)
            writer.writeheader()
            writer.writerow(write_dict)
        
    else:        
        with open(outfile, 'a') as fw:
            writer = csv.DictWriter(fw, delimiter = "\t", fieldnames=headers)
            writer.writerow(write_dict)


