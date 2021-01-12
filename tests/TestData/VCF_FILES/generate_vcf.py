#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import pandas
import contextlib
import unicodedata

def main(snakemake=None, cmd=None):
    parser = build_argument_parser()
    if snakemake is not None:
        cmd = get_cmd_from_snakemake(snakemake)
    args = parser.parse_args(cmd)
    run_vcf_generator(args)


def get_cmd_from_snakemake(snakemake):
    cmd = None
    return cmd

def build_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("mutation_file",
                        help="yaml describing mutations",
                        type=lambda f: validate_mutation_file(os.path.realpath(f)))
    parser.add_argument("--outfile", "-o")
    parser.add_argument('--mutation_ids', nargs="+")
    return parser


def validate_mutation_file(_file):
    _file = unicodedata.normalize("NFKC", _file).strip() # Just do not ask how this line ruinded an otherwise perfect day
    with open(_file, 'r') as mutfh:
        mutfile = yaml.load(mutfh, Loader=yaml.FullLoader)
    try:
        mutfile["reference"]
        mutfile["mutations"]
    except KeyError:
        raise ValueError("Mutation Yaml Does not conform to Spec")
    return mutfile

def run_vcf_generator(args):

    header = (
        f"##fileformat=VCFv4.2\n"
        f"##fileDate=20090805\n"
        f"##source=myImputationProgramV3.1\n"
        f'##reference={args.mutation_file["reference"]}'
    )
    struc_var_header=""

    data_header=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]

    data_header_str = f"\n#{chr(9).join(data_header)}"

    variation_data = gen_variation_data(args, data_header)

    with smart_open(args.outfile) as out_fh:
        print(f"{header}{struc_var_header}{data_header_str}", file=out_fh)
        variation_data.to_csv(out_fh, sep="\t",header=False, index=False)

def gen_variation_data(args, data_header):
    ref_name = args.mutation_file["reference"]
    mut_ref = args.mutation_file["mutations"]
    muts_selected = []
    if args.mutation_ids is None:
        muts_selected = list(mut_ref.keys())

    mut_df = pandas.DataFrame(dict(zip(data_header,
                                  [["."]*(len(muts_selected))]*len(data_header) )))
    mut_df["CHROM"] = ref_name
    mut_df["ID"] = muts_selected
    mut_df["POS"] = [mut_ref[key]["pos"] for key in muts_selected]
    mut_df["REF"] = [mut_ref[key]["ref"] for key in muts_selected]
    mut_df["ALT"] = [mut_ref[key]["alt"] for key in muts_selected]

    return mut_df

@contextlib.contextmanager
def smart_open(filename=None):
    """ Nice utillity open that can deal with stdout
    Thanks: Wolph
    https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
    """
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        main()
    else:
        main(snakemake=snakemake)
