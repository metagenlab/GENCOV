#!/usr/bin/env python


import pandas

if __name__ == '__main__':
    import argparse
    import sys
    import datetime
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", '--input', type=str, help="Input bed")
    parser.add_argument("-o", '--output', type=str, default="out.bed")

    args = parser.parse_args()

    df = pandas.read_csv(args.input, sep="\t", names=["chrom", "start", "end", "primer", "id", "strand", "seq"])



    '''

MN908947.3	25	50	SARS-CoV-2_1_LEFT	1	+	AACAAACCAACCAACTTTCGATCTC
MN908947.3	408	431	SARS-CoV-2_1_RIGHT	1	-	CTTCTACTAAGCCACAAGTGCCA
MN908947.3	324	344	SARS-CoV-2_2_LEFT	2	+	TTTACAGGTTCGCGACGTGC

=>

chrom	left_start	left_end	right_start	right_end
NC_045512.2	8	34	111	140
NC_045512.2	70	94	192	212
    '''


    # SARS-CoV-2_1_LEFT 
    df["location"] = df["primer"].apply(lambda x: x.split("_")[-1].lower())
    df["amplicon_id"] = df["primer"].apply(lambda x: int(x.split("_")[-2]))

    df_format = pandas.pivot_table(df, index=["chrom", "amplicon_id"], values=["start", "end"],columns=["location"]).reset_index().sort_values(["amplicon_id"])

    print(df_format.columns.tolist())

    df_format.columns = ['_'.join(i[::-1]) if i[1] != '' else i[0] for i in df_format.columns.tolist()]

    df_format[["chrom","left_start","left_end","right_start", "right_end"]].to_csv(args.output, sep="\t")