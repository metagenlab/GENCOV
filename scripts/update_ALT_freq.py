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
                

def update_alt_freq(in_fname, 
                    out_fname, 
                    min_freq,
                    ao_tag="AD", 
                    dp_tag="DP"):
    print("ao_tag", ao_tag)
    print("dp_tag", dp_tag)
    print("min_vf", min_freq)
    print("in_fname", in_fname)
    print("out_fname", out_fname)
    out_gz = out_fname.endswith(".gz")
    intermediate = re.sub("\.gz$","", out_fname)
    vcf_reader =  vcf.Reader(filename=in_fname)
    vcf_writer = vcf.Writer(open(intermediate, 'w'), vcf_reader)

    for call in vcf_reader:
        ref = call.REF
        alt = call.ALT[0]
        depth = int(call.INFO[dp_tag])
        ref_count, alt_count = call.samples[0][ao_tag]

        alt_fract = float(alt_count)/depth
        # skip call with frequency lower than min_freq
        if alt_fract < float(min_freq):
            continue
        # update AF field based on result
        call.INFO["AF"] = alt_fract

        vcf_writer.write_record(call)
    vcf_writer.flush()
    vcf_writer.close()
    if out_gz:
        print("out_gz", out_fname)
        bgzip_outname(intermediate, out_fname)


input_vcf_gz = snakemake.input[0]
output_vcf_gz =  snakemake.output[0]
ao_tag = snakemake.params["ao_tag"]
dp_tag = snakemake.params["dp_tag"]
min_freq = snakemake.params["min_freq"]
update_alt_freq(in_fname=input_vcf_gz,
                out_fname=output_vcf_gz,
                ao_tag=ao_tag,
                dp_tag=dp_tag,
                min_freq=min_freq)
