import os
import sys
import yaml
import argparse
from Bio import SeqIO
import csv
from collections import defaultdict

sample_names = ["denv1_test", "denv2_test", "denv3_test", "denv4_test", "low_coverage"]
correct_serotype = {"denv1_test": ["DENV1"], "denv2_test":["DENV2"],  "denv3_test":["DENV3"],  "denv4_test":["DENV4", "DENV4_sylvatic"]}
correct_top_call = {"denv1_test": "DENV1", "denv2_test":"DENV2",  "denv3_test":"DENV3",  "denv4_test": "DENV4", "low_coverage":"NA"}
correct_sero_to_seq = {'DENV1': 'denv1_test', 'DENV2': 'denv2_test', 'DENV3': 'denv3_test', 'DENV4': 'denv4_test', 'DENV4_sylvatic': 'denv4_test'}
serotypes = ["DENV1", "DENV2", "DENV3", "DENV4", "DENV4_sylvatic"]

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def check_alignments_present(outdir):

    for sero in serotypes:
        assert os.path.exists(os.path.join(outdir, "alignments", f"{sero}.trim.aln"))
        assert os.path.exists(os.path.join(outdir, "alignments", f"{sero}.untrim.aln"))

    return

def check_alignments_correct(outdir):

    alignment_dir = os.path.join(outdir, "alignments")

    for sero in serotypes:
        correct_sequence = f'{correct_sero_to_seq[sero]}.{sero}'
        for seq in SeqIO.parse(os.path.join(alignment_dir, f"{sero}.trim.aln"), 'fasta'):
            assert seq.id == correct_sequence
        for seq in SeqIO.parse(os.path.join(alignment_dir, f"{sero}.untrim.aln"), 'fasta'):
            assert seq.id == correct_sequence

    return

def check_files_present(outdir, dir_name, file_end, file_end2=False, other_file=False):

    files_present = []

    for sample in sample_names:
        if sample in correct_serotype:
            for sero in correct_serotype[sample]:
                rel_file = os.path.join(outdir, dir_name, f'{sample}.{sero}.{file_end}')
                files_present.append(rel_file)
                print(rel_file)
                assert os.path.exists(rel_file)
                assert os.stat(rel_file).st_size > 0

                if file_end2:
                    rel_file2 = os.path.join(outdir, dir_name, f'{sample}.{sero}.{file_end2}')
                    files_present.append(rel_file2)
                    assert os.path.exists(rel_file2)
                    assert os.stat(rel_file2).st_size > 0

    if other_file:
        other_file_path = os.path.join(outdir, dir_name, other_file)
        files_present.append(other_file_path)
        assert os.path.exists(other_file_path)
        assert os.stat(other_file_path).st_size > 0

    for f in os.listdir(os.path.join(outdir, dir_name)):
        try:
            assert os.path.join(outdir, dir_name, f) in files_present
        except AssertionError:
            print(os.path.join(outdir, dir_name, f))
            print(files_present)

    return


def check_bams_present(outdir):

    directory = "bam_files"
    file_end = "sort.bam"
    file_end2 = "sort.bam.bai"

    check_files_present(outdir, directory, file_end, file_end2=file_end2)

    return


def check_consensus_present(outdir, depth):

    directory = "consensus_sequences"
    file_end = f'{depth}.cons.fa'    
    
    check_files_present(outdir, directory, file_end)

    return


def check_consensus_correct(outdir, depth):

    dir_path = os.path.join(outdir, "consensus_sequences")

    for sample in sample_names:
        if sample in correct_serotype:
            for sero in correct_serotype[sample]:
                correct_sequence = f'{sample}.{sero}'
                for seq in SeqIO.parse(os.path.join(dir_path, f"{sample}.{sero}.{depth}.cons.fa"), 'fasta'):
                    assert seq.id == correct_sequence
    
    return



def check_depth_present(outdir):

    directory = "depth"
    file_end = "depth.txt"

    check_files_present(outdir, directory, file_end)

    return


def check_depth_correct(outdir):

    correct_lengths ={'denv1_test.DENV1': 10735, 'denv2_test.DENV2': 10723, 'denv3_test.DENV3': 10707, 'denv4_test.DENV4': 10649, 'denv4_test.DENV4_sylvatic': 10667}

    for sample in sample_names:
        if sample in correct_serotype:
            for sero in correct_serotype[sample]:
                file = os.path.join(outdir, "depth", f'{sample}.{sero}.depth.txt')
                with open(file) as f:
                    data = csv.DictReader(f, delimiter="\t")
                    assert data.fieldnames == ["position", "depth"]

                with open(file) as f:
                    assert sum(1 for row in f) -1 == correct_lengths[f'{sample}.{sero}']

    return

def check_variants_present(outdir, depth):

    directory = "variants"
    file_end = f"{depth}.variants_frequency.tsv"

    summary_file = "variants_summary.tsv"

    check_files_present(outdir, directory, file_end, other_file=summary_file)

    return


def check_variants_correct(outdir, depth):

    correct_headers = ['position', 'reference_base', 'alternative_base', 'reference_depth', 'reference_depth_reverse', 'reference_quality', 'alternate_depth', 'alternate_depth_reverse', 'alternative_quality', 'alternative_frequency', 'total_depth', 'p_value_fisher', 'pass', 'gff_feature', 'reference_codon', 'reference_amino_acid', 'alternative_codon', 'alternative_amino_acid']
    digit_cols = ['position', 'reference_depth', 'reference_depth_reverse', 'reference_quality', 'alternate_depth', 'alternate_depth_reverse', 'alternative_quality', 'alternative_frequency', 'total_depth', 'p_value_fisher']

    num_variants = {}

    for sample in sample_names:
        if sample in correct_serotype:
            for sero in correct_serotype[sample]:
                file = os.path.join(outdir, "variants", f'{sample}.{sero}.{depth}.variants_frequency.tsv')
                with open(file) as f:
                    data = csv.DictReader(f, delimiter="\t")
                    assert data.fieldnames == correct_headers
                
                with open(file) as f:
                    num_variants[f'{sample}.{sero}'] = sum(1 for row in f) -1 
 
                with open(file) as f:
                    data = csv.DictReader(f, delimiter="\t")
                    for l in data:
                        for col in digit_cols:
                            assert isfloat(l[col])

    return num_variants


def check_variant_summary(outdir, num_variants):

    correct_headers = ["sample_id", "serotype", "variant_count"]
    
    with open(os.path.join(outdir, "variants", "variants_summary.tsv")) as f:
        data = csv.DictReader(f, delimiter="\t")
        assert data.fieldnames == correct_headers
        for l in data:
            assert l['sample_id'] in sample_names
            assert l['serotype'] in correct_serotype[l['sample_id']]
            assert int(l['variant_count']) == num_variants[f'{l["sample_id"]}.{l["serotype"]}']

    return

def check_summaries_present(outdir):

    
    for sum in ["summary_all_samples.tsv", "top_virus_all_samples.tsv", "virus_calls.tsv"]:

        assert os.path.exists(os.path.join(outdir, sum))
        assert os.stat(os.path.join(outdir, sum)).st_size > 0

    return


def check_all_summary_correct(outdir):

    sum_file = os.path.join(outdir, "summary_all_samples.tsv")
    headers = ['sample_id', 'consensus_sequence_file', 'depth', 'serotype_called', 'reference_sequence_name', 'reference_sequence_length', 'number_aligned_bases', 'coverage_untrimmed', 'coverage_trimmed']
    digit_cols = ['depth', 'reference_sequence_length', 'number_aligned_bases', 'coverage_untrimmed', 'coverage_trimmed']

    num_lines = {}

    with open(sum_file) as f:
        data = csv.DictReader(f, delimiter="\t")
        assert data.fieldnames == headers
        for l in data:

            if l["sample_id"] not in num_lines:
                num_lines[l['sample_id']] = 1
            else:
                num_lines[l['sample_id']] += 1

            for col in digit_cols:
                assert isfloat(l[col]) or "NA"

    assert len(num_lines) == len(sample_names)
    
    for k,v in num_lines.items():
        assert v == 6

    return


def check_top_calls_correct(outdir):

    headers = ['sample_id', 'consensus_sequence_file', 'depth', 'serotype_called', 'reference_sequence_name', 'reference_sequence_length', 'number_aligned_bases', 'coverage_untrimmed', 'coverage_trimmed']
    num_lines = {}
    top_call = {}
    with open(os.path.join(outdir, "top_virus_all_samples.tsv")) as f:
        data = csv.DictReader(f, delimiter="\t")
        assert data.fieldnames == headers
        for l in data:
            if l["sample_id"] not in num_lines:
                num_lines[l['sample_id']] = 1
            else:
                num_lines[l['sample_id']] += 1

            top_call[l['sample_id']] = l['serotype_called']

    assert len(num_lines) == len(sample_names)
    
    for k,v in num_lines.items():
        assert v == 1
    for k,v in top_call.items():
        assert v == correct_top_call[k]

    return

def check_virus_calls_correct(outdir):
    
    headers = ['sample_id', 'consensus_sequence_file', 'depth', 'serotype_called', 'reference_sequence_name', 'reference_sequence_length', 'number_aligned_bases', 'coverage_untrimmed', 'coverage_trimmed']

    correct_calls = defaultdict(list)
    with open(os.path.join(outdir, "virus_calls.tsv")) as f:
        data = csv.DictReader(f, delimiter="\t")
        assert data.fieldnames == headers
        for l in data:
            correct_calls[l['sample_id']].append(l['serotype_called'])

    assert correct_calls == correct_serotype

    return


def run_tests():

    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir")
    parser.add_argument("--depth")

    args = parser.parse_args()
    outdir = args.outdir
    depth = args.depth

    check_alignments_present(outdir)
    check_alignments_correct(outdir)

    check_bams_present(outdir)

    check_consensus_present(outdir, depth)
    check_consensus_correct(outdir, depth)

    check_depth_present(outdir)
    check_depth_correct(outdir)

    check_variants_present(outdir, depth)
    num_variants = check_variants_correct(outdir, depth)
    check_variant_summary(outdir, num_variants)

    check_summaries_present(outdir)
    check_all_summary_correct(outdir)
    check_top_calls_correct(outdir)
    check_virus_calls_correct(outdir)


if __name__=="__main__":
    run_tests()