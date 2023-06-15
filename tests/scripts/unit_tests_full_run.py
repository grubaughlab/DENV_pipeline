import os
import sys
import yaml
import argparse
from Bio import SeqIO


sample_names = ["denv1_test", "denv2_test", "denv3_test", "denv4_test", "low_coverage"]
correct_serotype = {"denv1_test": ["DENV1"], "denv2_test":["DENV2"],  "denv3_test":["DENV3"],  "denv4_test":["DENV4", "DENV4_sylvatic"] }
correct_sero_to_seq = {'DENV1': 'denv1_test', 'DENV2': 'denv2_test', 'DENV3': 'denv3_test', 'DENV4': 'denv4_test', 'DENV4_sylvatic': 'denv4_test'}
serotypes = ["DENV1", "DENV2", "DENV3", "DENV4", "DENV4_sylvatic"]

def check_alignments_present(outdir):

    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV1.trim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV1.untrim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV2.trim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV2.untrim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV3.trim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV3.untrim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV4.trim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV4.untrim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV4_sylvatic.trim.aln"))
    assert os.path.exists(os.path.join(outdir, "results", "alignments", "DENV4_sylvatic.untrim.aln"))

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


def check_bams_present(outdir):

    for sample in sample_names:
        if sample in correct_serotype:
            for sero in correct_serotype[sample]:

                bam_file = os.path.join(outdir, "results", "bam_files", f"{sample}.{sero}.sort.bam")
                bai_file = os.path.join(outdir, "results", "bam_files", f"{sample}.{sero}.sort.bam.bai")

                assert os.path.exists(bam_file)
                assert os.path.exists(bai_file)

                assert os.stat(bam_file).st_size > 0
                assert os.stat(bai_file).st_size > 0

    return


def check_consensus_present(outdir):

    return


def check_consensus_correct(outdir):


    return

def check_depth_present(outdir):

    return


def check_depth_correct(outdir):


    return

def check_variants_present(outdir):

    return


def check_variants_correct(outdir):


    return

def check_summaries_present(outdir):

    return


def check_summaries_correct(outdir):


    return


def run_tests():

    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir")

    args = parser.parse_args()
    outdir = args.outdir

    check_alignments_present(outdir)
    check_alignments_correct(outdir)

    check_bams_present(outdir)

    check_consensus_present(outdir)
    check_consensus_correct(outdir)

    check_depth_present(outdir)
    check_depth_correct(outdir)

    check_variants_present(outdir)
    check_variants_correct(outdir)

    check_summaries_present(outdir)
    check_summaries_correct(outdir)


if __name__=="__main__":
    run_tests()