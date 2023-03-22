import os
import csv
import argparse
from collections import defaultdict
import shutil

def summarise_files(config, serotype_calls):

    # parser = argparse.ArgumentParser()
    # parser.add_argument("--serotype-calls", dest="serotype_calls")
    # parser.add_argument("--config")

    # args = parser.parse_args()
    # config = args.config

    min_coverage = []
    serotypes = defaultdict(list)
    with open(serotype_calls) as f:
        data = csv.DictReader(f, delimiter="\t")
        for l in data:
            if int(l['CoverageUntrimmed']) >= 50:
                min_coverage.append(l)
                serotypes[l['SampleID']].append(l['Serotype'])

    sort_variant_files(config, serotypes)
    get_right_serotype_files(config, serotypes)
    
    if config["temp"]:
        with open(os.path.join(config["tempdir"], "DENV.serotype.calls.mincov50.final.tsv"), 'w') as fw:
            headers = ["sample_id", "consensus_sequence_file", "depth", "serotype", "reference_serotype_name", "reference_sequence_length", "number_aligned_bases", "coverage_untrimmed", "coverage_trimmed"]
            writer = csv.DictWriter(fw, delimiter="\t")
            writer.writeheader(fieldnames=headers)
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
        if file.endswith("_variants.tsv"):
            count = 0
            new_file_name = file.replace("_variants.tsv", "variants_frequency.tsv")
            new_file = os.path.join(config['outdir'], new_file_name)
            with open(new_file, 'w') as fw:
                writer = csv.DictWriter(fw, delimiter="\t")
                writer.writeheader(fieldnames=list(old_to_new.values()))
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
    with open(new_depth_file) as fw:
        writer = csv.DictWriter(fw, delimiter="\t")
        writer.writeheader(fieldnames=headers)
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
    for sample, serotype_lst in serotypes.items():
        for serotype in serotype_lst:
            bam_file = f'{sample}.{serotype}.sort.bam'
            consensus_file = f'{sample}.{serotype}.cons.fa'
            depth = f'{sample}.{serotype}.depth.txt'
            variant_frequency = f'{sample}.{serotype}.{config["depth"]}_variants_frequency.tsv'

            bam_files.add(bam_file)
            consensus.add(consensus_file)
            depths.add(depth)
            variant_frequencies.add(variant_frequency)

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

    for var_freq in variant_frequencies:
        source = os.path.join(config['outdir'], var_freq)
        dest = os.path.join(config["outdir"], "results", "variants")
        shutil.move(source, dest)


