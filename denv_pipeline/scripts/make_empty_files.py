
import os
import csv

import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--depth")
    parser.add_argument("--tempdir")
    parser.add_argument("--sample-name", dest="sample_name")
    parser.add_argument("--virus-type", dest="virus_type")

    args = parser.parse_args()
    
    tempdir = args.tempdir
    depth = args.depth
    sample_name = args.sample_name
    virus_type = args.virus_type

    do_virus_summary_file(tempdir, sample_name, depth, virus_type)
    do_variant_file(tempdir, sample_name, depth, virus_type)


def do_virus_summary_file(tempdir, sample_name, depth, virus_type):
    
    headers = ["sample_id","consensus_sequence_file","depth","serotype_called","reference_sequence_name","reference_sequence_length","number_aligned_bases","coverage_untrimmed","coverage_trimmed"]

    outfile = os.path.join(tempdir, f'{sample_name}_all_virustype_info.txt')

    write_dict = {}
    write_dict["sample_id"] = sample_name
    write_dict["consensus_sequence_file"] = "NA"
    write_dict['depth'] = depth

    write_dict["serotype_called"] = "NA"
    write_dict["reference_sequence_name"] = virus_type
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


def do_variant_file(tempdir, sample_name, depth, virus_type):

    outfile = os.path.join(tempdir, f"{sample_name}.{virus_type}.{depth}.variants.tsv")

    headers = ["REGION","POS","REF","ALT","REF_DP","REF_RV","REF_QUAL","ALT_DP","ALT_RV","ALT_QUAL","ALT_FREQ","TOTAL_DP","PVAL","PASS","GFF_FEATURE","REF_CODON","REF_AA","ALT_CODON","ALT_AA","POS_AA"]

    with open(outfile, 'w') as fw:
        writer = csv.DictWriter(fw, delimiter = "\t", fieldnames=headers)
        writer.writeheader()

if __name__=="__main__":
    main()

    