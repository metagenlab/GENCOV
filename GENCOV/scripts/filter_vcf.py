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


                

def update_alt_freq(in_fname, 
                    out_fname, 
                    min_freq,
                    modif_pos_out,
                    keep_variant_list):

    keep_variant_list = [int(i) for i in keep_variant_list]

    vcf_reader =  vcf.Reader(filename=in_fname)
    vcf_writer = vcf.Writer(open(out_fname, 'w'), vcf_reader)
    pos_out = open(modif_pos_out, 'w')
    pos_out.write("post\tref\talt\talt_freq\tEPPR\n")
    
    for call in vcf_reader:
        alt_fract = float(call.INFO["AF"][0])
        # skip call with frequency lower than min_freq
        if alt_fract <= float(min_freq):
            # remove variant below min_freq, except if among known problematic positions
            if int(call.POS) not in keep_variant_list:
                continue
            else:
                print(f"KEEPING:\t{call.REF}{call.POS}{call.ALT}\tfreq: {alt_fract}")
                pos_out.write(f'{call.POS}\t{call.REF}\t{call.ALT[0]}\t{alt_fract}\t{call.INFO["EPPR"]}\n')

        vcf_writer.write_record(call)
    vcf_writer.flush()
    vcf_writer.close()
    pos_out.close()

input_vcf_gz = snakemake.input[0]

output_vcf_gz =  snakemake.output[0]
modified_pos = snakemake.output[1]

min_freq = snakemake.params["frac_filter"]
position_list = snakemake.params["keep_positions"]

print("--- FILTER VCF ---")
print("KEEP", position_list)

update_alt_freq(in_fname=input_vcf_gz,
                out_fname=output_vcf_gz,
                min_freq=min_freq,
                modif_pos_out=modified_pos,
                keep_variant_list=position_list)
