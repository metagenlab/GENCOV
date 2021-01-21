
import pandas 
import os
from Bio.SeqUtils import seq1
import re 
variant_tables = snakemake.input 
import numpy as np

sample2mutations = {}

aa_map = {'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', 'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y', 'Ter': '*'}

def multiple_replace(dict, text):
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text) 


df_list = []
for variant_table in variant_tables:
    sample = os.path.basename(variant_table).split(".")[0]
    df = pandas.read_csv(variant_table, sep='\t')
    print(df.columns)
    df.columns = ["CHROM", "POS","REF","ALT","EFFECT","GENE","Protein_Change","AA_POS", "Depth", "RefNum", "AltNum"]

    df["sample"] = [sample] * len(df.index)

    updated_change_aa = []
    updated_change_aa_with_prot = []
    updated_change_nucl = []
    alt_percent = []
    for n,row in df.iterrows():
        aa_change = row["Protein_Change"]
        print(aa_change)
        # replace three letter code by single letter code
        # p.Gln62_Phe63delinsHis => p.Q62_F63delinsH
        # p.His83fs
        if not pandas.isna(aa_change):
            single_letter = re.sub("p\.", "", multiple_replace(aa_map, aa_change))
            updated_change_aa.append(single_letter)
            updated_change_aa_with_prot.append(f'{single_letter} ({row["GENE"]})')
        else:
            updated_change_aa.append('n/a')
        # nucl change
        updated_change_nucl.append(f'{row["REF"]}{row["POS"]}{row["ALT"]}')

        alt_percent.append(round((float(row["AltNum"])/float(row["Depth"]))*100,2))
    df["AA_change"] = updated_change_aa
    df["Nucl_change"] = updated_change_nucl
    df["Var_percent"] = alt_percent

    sample2mutations[sample] = {}
    sample2mutations[sample]["aa_changes"] = ';'.join([i for i in updated_change_aa_with_prot])
    sample2mutations[sample]["nucl_changes"] = ';'.join(updated_change_nucl)

    df_list.append(df)

combined_df = pandas.concat(df_list)
df_format = combined_df[ ["sample","CHROM", "POS","REF","ALT","EFFECT","GENE", "AA_POS", "AA_change", "Nucl_change", "Depth", "RefNum", "AltNum","Var_percent"]]

df_format.to_csv(snakemake.output[0], sep='\t', index=False)



with open(snakemake.output[1], 'w') as f:
    for sample in set(combined_df["sample"]):
        f.write(f"{sample}\t{sample2mutations[sample]}\n")
