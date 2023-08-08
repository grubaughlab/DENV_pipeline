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
        possible_tops = []
        with open(file) as f:
            data = csv.DictReader(f, delimiter="\t")
            for l in data:
                possible_tops.append(l)
                all_lines.append(l)

                if float(l['coverage_untrimmed']) >= 50:
                    serotype_call.append(l)
                    serotypes[l['sample_id']].append(l['serotype_called'])
            
        if len(possible_tops) > 0:
            top = sorted(possible_tops, key=lambda x:float(x['coverage_untrimmed']), reverse=True)[0]
            top_calls.append(top)
        elif len(possible_tops) == 1:
            top_calls.append(possible_tops[0])
        else:
            pass

    headers = ["sample_id","consensus_sequence_file","depth","serotype_called","reference_sequence_name","reference_sequence_length","number_aligned_bases","coverage_untrimmed","coverage_trimmed"]
    
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
    get_right_serotype_files(config, serotypes)
    make_alignments(config, serotypes)

    return


def sort_variant_files(config, serotypes):

    old_to_new = {'POS': 'position', 'REF': 'reference_base', 'ALT': 'alternative_base', 'REF_DP': 'reference_depth', 'REF_RV': 'reference_depth_reverse', 'REF_QUAL': 'reference_quality', 'ALT_DP': 'alternate_depth', 'ALT_RV': 'alternate_depth_reverse', 'ALT_QUAL': 'alternative_quality', 'ALT_FREQ': 'alternative_frequency', 'TOTAL_DP': 'total_depth', 'PVAL': 'p_value_fisher', 'PASS': 'pass', 'GFF_FEATURE': 'gff_feature', 'REF_CODON': 'reference_codon', 'REF_AA': 'reference_amino_acid', 'ALT_CODON': 'alternative_codon', 'ALT_AA': 'alternative_amino_acid'}

    summary_file = open(os.path.join(config["outdir"], "results", "variants", "variants_summary.tsv"), 'w')
    summary_file.write("sample_id\tserotype\tvariant_count\n")

    full_virus_type_list = config["virus_type_list"]
    depth = config["depth"]
    for sample, serotype_lst in serotypes.items():
        for option in full_virus_type_list:
            input_file = os.path.join(config["tempdir"], f"{sample}.{option}.{depth}.variants.tsv")
            output_file = os.path.join(config["tempdir"], f"{sample}.{option}.{depth}.variants_frequency.tsv")

            count = 0
            with open(output_file, 'w') as fw:
                headers = list(old_to_new.values())
                writer = csv.DictWriter(fw, headers, delimiter="\t")
                writer.writeheader()

                with open(input_file) as f:
                    data = csv.DictReader(f, delimiter="\t")
                    for l in data: 
                        write_dict = {}
                        if l['PASS'] == "TRUE" and float(l['ALT_FREQ']) < 0.8 and float(l['ALT_FREQ']) > 0.2:
                            for old,new in old_to_new.items():
                                write_dict[new] = l[old]
                            writer.writerow(write_dict)
                            count += 1

            if option in serotype_lst:
                summary_file.write(f"{sample}\t{option}\t{count}\n")

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
    depth = config["depth"]
    full_virus_type_list = config["virus_type_list"]
    
    for sample, serotype_lst in serotypes.items():
        for option in serotype_lst:
            
            bam_files.add(f'{sample}.{option}.sort.bam')
            bam_indices.add(f'{sample}.{option}.sort.bam.bai')
            consensus.add(f'{sample}.{option}.{depth}.cons.fa')
            depths.add(f'{sample}.{option}.depth.txt')
            variant_frequencies.add(f'{sample}.{option}.{depth}.variants_frequency.tsv')


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


def make_alignments(config, serotypes):

    alignment_dir = os.path.join(config["outdir"], "results", "alignments")
    tempdir = config["tempdir"]
    depth = config["depth"]
    full_virus_type_list = config["virus_type_list"]
    
    trimmed_aln_dict = defaultdict(list)
    untrimmed_aln_dict = defaultdict(list)
    for sample, serotype_lst in serotypes.items():
        for i in serotype_lst:
            trimmed = f'{sample}.{i}.{depth}.out.trim.aln'
            untrimmed = f'{sample}.{i}.{depth}.out.aln'
            
            trimmed_aln_dict[i].append(trimmed)
            untrimmed_aln_dict[i].append(untrimmed)

    for virus_type, alns in trimmed_aln_dict.items():
        new_file = f'{alignment_dir}/{virus_type}.trim.aln'
        cat_string = " ".join([f"{tempdir}/{aln}" for aln in alns])
        os.system(f"cat {cat_string} >> {new_file}")

    for virus_type, alns in untrimmed_aln_dict.items():
        new_file = f'{alignment_dir}/{virus_type}.untrim.aln'
        cat_string = " ".join([f"{tempdir}/{aln}" for aln in alns])
        os.system(f"cat {cat_string} >> {new_file}")
