#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF-1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import argparse
import subprocess
import re
import sys
import gzip
import vcf 


def parse_args(CMD=None):
    parser = argparse.ArgumentParser(prog="rename_in_gff3.py", description="changes genotype in VCFs", )
    parser.add_argument('vcf', metavar="FILE", help="vcf file", type=str)
    parser.add_argument('--ao', metavar="STR", help="tag for read count supporting the respective variant (default: AO)", type=str, default="AD")
    parser.add_argument('--dp', metavar="STR", help="tag for total read count at the repsective position (default: DP)", type=str, default="DP")
    parser.add_argument('-o', help="output file (will be overwritten!)", type=str, required=True)
    parser.add_argument('--vf', metavar="FLOAT", help="minimal variant fraction to set a homogeneous genotype (default: 0.9)", type=float, default=0.9)
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    return parser.parse_args(CMD)

def get_cmd_from_snake(snakemake):
    cmd = ["-o", snakemake.output[0],
           "--vf", snakemake.params.frac,
          snakemake.input[0]]
    cmd = [str(arg) for arg in cmd]
    return cmd
                
def process(in_fname, out_fname, min_vf, ao_tag="ao", dp_tag="AO"):
    print("ao_tag", ao_tag)
    print("dp_tag", dp_tag)
    print("min_vf", min_vf)
    print("in_fname", in_fname)
    print("out_fname", out_fname)
    #sanity checks
    if min_vf <= 0.5:
        sys.exit("error: min_vf has to be greater than 0.5")
    out_gz = out_fname.endswith(".gz")
    intermediate = re.sub("\.gz$","", out_fname)

    vcf_reader =  vcf.Reader(filename=in_fname)
    vcf_writer = vcf.Writer(open(intermediate, 'w'), vcf_reader)
    keep=[]
    for call in vcf_reader:
        print(call)
        ref = call.REF
        alt = call.ALT[0]
        depth = int(call.INFO[dp_tag])
        ref_count, alt_count = call.samples[0][ao_tag]
        alt_fract = float(alt_count)/depth
        print(alt_fract)
        if alt_fract < min_vf:
            print(f"{alt_fract} smaller than {min_vf}")
            continue
        print(f"Writing {ref}{call.POS}{alt}")
        vcf_writer.write_record(call)
    vcf_writer.flush()
    vcf_writer.close()
    if out_gz:
        print("out_gz", out_fname)
        bgzip_outname(intermediate, out_fname)
        
def bgzip_outname(_file, outfile=None):
    if outfile is not None:
        with open(outfile, "wb") as out_fh:
            with subprocess.Popen(["bgzip", _file, "-c"], 
                                   stdout=out_fh) as bg_proc:
                ret = bg_proc.wait() 
        os.remove(_file)  
    else:
        with subprocess.Popen(["bgzip", _file], 
                               stderr=subprocess.PIPE) as bg_proc:
            out, err = bg_proc.communicate()
            ret = bg_proc.wait() 
                
def main(CMD=None):
        args = parse_args(CMD)
        process(args.vcf, args.o, args.vf, args.ao, args.dp)

if __name__ == "__main__":
        CMD = None
        if "snakemake" in globals(): 
            CMD = get_cmd_from_snake(snakemake)
        print("CMD", CMD)
        main(CMD=CMD)
