import os
import sys
import yaml
import argparse
from Bio import SeqIO


sample_names = ["denv1_test", "denv2_test", "denv3_test", "denv4_test", "low_coverage"]
correct_serotype = {"denv1_test": ["DENV1"], "denv2_test":["DENV2"],  "denv3_test":["DENV3"],  "denv4_test":["DENV4", "DENV4_sylvatic"] }

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





    return


def check_bams_present(outdir):



    return


def check_bams_correct(outdir):


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
    check_bams_correct(outdir)

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