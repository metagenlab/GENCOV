#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import random

def main(snakemake=None, cmd=None):
    parser = build_argument_parser()
    if snakemake is not None:
        cmd = get_cmd_from_snakemake(snakemake)
    args, unknown_args = parser.parse_known_args(cmd)
    valid_setup = validate_setup(args,unknown_args)
    finalize_setup(valid_setup)

def get_cmd_from_snakemake(snakemake):
    cmd = None
    return cmd


def build_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("guide_file",
                        help="yaml descibing what needs to be generated",
                        type=lambda f: os.path.realpath(f))
    parser.add_argument("--outfolder", "-o")
    return parser


def validate_setup(args, unknown_args):
    with open(args.guide_file, 'r') as fh_guide:
        test_yaml_config = yaml.load(fh_guide, Loader=yaml.FullLoader)
        #import pdb; pdb.set_trace()
    test_yaml_config = secondary_argparse(test_yaml_config, unknown_args)
    return test_yaml_config


def secondary_argparse(test_setup, cmd):
    setup = test_setup['test']
    tbd_fields = [ key for key in setup if setup[key] in ["tbd", ""]]
    custom_parse_items = []
    if "reference" in tbd_fields:
        custom_parse_items.append("reference")
    parser = argparse.ArgumentParser()
    parse_ref = {"reference": [["--reference"],
                               {"required":True,
                                "type":lambda f: os.path.realpath(f)}]}

    for parse_item in custom_parse_items:
        parser.add_argument(*parse_ref[parse_item][0],**parse_ref[parse_item][1])
    args = parser.parse_args(cmd)
    setup.update(vars(args))
    print(test_setup)
    return test_setup

def finalize_setup(valid_setup):
    setup = add_necessary_simulation_runs(valid_setup)
    prepare_testfolder(setup)

def add_necessary_simulation_runs(valid_setup):
    test_setup = valid_setup["test"]
    seed = test_setup["seed"]
    replicate_setup = {}
    tmp_seed = seed
    random.seed(a=tmp_seed)
    for rep_num in range(test_setup["replicates"]):
        rep_id = f"rep_{rep_num + 1}"
        replicate_setup[rep_id] = {}
        tmp_seed = random.randint(1, sys.maxsize)
        random.seed(a=tmp_seed)
        replicate_setup[rep_id]["seed"] = tmp_seed

    print(replicate_setup)



def prepare_testfolder(setup):
    pass


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        main()
    else:
        main(snakemake=snakemake)
