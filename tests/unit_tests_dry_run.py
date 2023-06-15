import os
import sys
import yaml
import argparse

def file_structure_exists(outdir):

    assert os.path.exists(os.path.join(outdir, "downloads")) 
    assert os.path.exists(os.path.join(outdir, "log_files")) 
    assert os.path.exists(os.path.join(outdir, "test_temp"))
    assert os.path.exists(os.path.join(outdir, "results")) 
    assert os.path.exists(os.path.join(outdir, "output_config.yml"))

    assert os.path.exists(os.path.join(outdir, "results", "alignments"))
    assert os.path.exists(os.path.join(outdir, "results", "bam_files"))
    assert os.path.exists(os.path.join(outdir, "results", "consensus_sequences"))
    assert os.path.exists(os.path.join(outdir, "results", "depth"))
    assert os.path.exists(os.path.join(outdir, "results", "variants"))


def args_correctly_added_from_config_file(in_config, out_config):

    with open(in_config,"r") as f:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
    
    with open(out_config, 'r') as f:
        output_config = yaml.load(f, Loader=yaml.FullLoader)
        
    assert input_config == output_config

def run_dryrun_tests():

    parser = argparse.ArgumentParser()

    parser.add_argument("--in-config", dest="in_config")
    parser.add_argument("--out-config", dest="out_config")
    parser.add_argument("--outdir")

    args = parser.parse_args()

    in_config = args.in_config
    out_config = args.out_config
    outdir = args.outdir

    file_structure_exists(outdir)
    args_correctly_added_from_config_file(in_config, out_config)


if __name__=="__main__":
    run_dryrun_tests()