# GENCOV: pipeline adapted from nCoV_minipipe (aka CovPipe)

[[_TOC_]]

## Introduction

CovPipe is a pipeline to generate consensus sequences from NGS reads based on a reference sequence.
The pipeline is tailored to be used for SARS-CoV-2 data, but may be used for other viruses.

Genomic variants of your NGS data in comparison to a reference will be determined.
These variants will be included into the reference and form the consensus sequences.
See below for further details on the determined set of consensus sequences.

## Setup

Start with fresh clone of the pipeline.
Please use one of our releases to achieve maximal performance.

https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe/-/releases

Alternatively, one can use a clone of the repository.

```bash
cd folder/you/want/to/setup/the/tool/in

git clone https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe.git
```


## Usage

There are two ways of executing the pipeline. 
While the use of a wrapper is recommended, a manual option is provided, if there should be issues with running the wrapper on your system.

### Using the wrapper

The wrapper has a few dependencies.

| Dependencies | version |
| ------------ | ------- |
| [python](https://www.python.org/downloads/) | 3.6.0+ |
| [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) | 4.5+ |
| [strictyaml](https://pypi.org/project/strictyaml/) | 1.0.6 |

You don't need to worry too much about the dependencies as all of them can be installed quite quickly using the "Quick setup with conda". This needs [conda](https://docs.conda.io/en/latest/miniconda.html) installed, though.


#### Quick setup with conda

There is a quick way to get the pipeline running.

Simply run (this needs conda installed on your system) 

```bash
bash ncov_minipipe.conda.setup/setup_env.sh [optional: path/where/to/put/covpipe_environment]
```

If no path is provided as parameter the script defaults to put the conda environment into the
current working directory. 
Optionally other locations can be used.

```bash
# if no prefix was provided should yield
ls .

covpipe_environment [other files] [and folders]
```

You then should be able to run the pipeline. 

#### Run pipeline in conda environment

Environment activation.
This step has to be executed only once for each analyses session.

```bash
source activate covpipe_environment/   # Do not forget the trailing slash
```

This results in a modified prompt of your console by adding "(tool-env: ncov_minipipe)".
This environment is used to control environment tools for the pipeline.

simply run:

```bash
ncov_minipipe --help
```

example run of the pipeline:

```bash
(tool-env: ncov_minipipe)$ ncov_minipipe \
--outputfolder ... \
--reference .../NC_045512.2.fasta \
--kraken .../GRCh38.p13_GBcovid19-2020-05-22 \
--taxid 2697049 \
--amplicon \
--primerfasta .../cleanplex.primer.cutadapt.v2.fasta \
--adapter adapter.fasta \
...folder/*fastq.gz
```

Optional argumets can be used to modify the behaviour of the pipeline.

- --run_id makes an identifier appear in the QC report, which may facilitate overview, if multiple data sets are processed
- --low\_af, --good\_af, --low\_qual and min\_depth are used to influence the variant filter steps

### Manual Way

The manual way requires that you collect all sample IDs and the corresponding data as well as configure the pipeline. 
Although this sounds difficult, it's done quite fast and template files are provided.

This way also works with python 2.7 on your system.

There is no need to activate any conda environments prior to the start of the pipeline to use this method.

1. create a sample yaml according to the template (template_sample.config.yaml). Make sure to use full paths.
2. modify the ncov_minipipe.config to contain paths to tools, sample config and adapter.fasta
    - it's best, if your fasta files do contain one of the suffixes .fasta, .fn, .fna or .fa
    - [optional: amplification primer removal] make sure the file used for amplification primer removal contains the suffix .fasta
3. run pipeline with a command line like
    - `$ snakemake -s ncov_minipipe.snake --configfile ncov_minipipe.config --cores 60 --use-conda`
    - make sure to run at least snakemake 5, which is available in the covpipe environment
    - Requires conda to run (obtainable here: [miniconda](https://docs.conda.io/en/latest/miniconda.html))
    - Needs network connection (possible proxy íssues) to download packages from conda repositories, e.g. bioconda and conda-forge
    - Please note that the first run of the pipeline may take considerably longer as the package install needs some time, but the packages are reused in later runs

## Files you need

- seq_ID_[12].fastq.gz : file(s) of reads (paired end)
- genome.fasta : genome reference file
- adapter.fasta : list of Illumina sequencing adapters that may show up at the 3' end of the reads
- ncov_minipipe.config : file containing all settings to the pipeline (template provided)
- sample_config.yaml : file containing all read file locations with full paths (template provided)
- optional:
- primer_sequences.fasta : list of primers that have been used in amplification. The sequences are to be found at the 5' end of the reads
- kraken2 database for read filtering. In case of SARS-CoV-2 data you can download a database to remove especially human reads from https://zenodo.org/record/3854856
  - please use taxID 2697049 together with the database

## Outputs

- results will be put to Analyses folder
- logs into logs folder
- the final report will locate to Documentation

## Pipeline steps

### Step 1: trimming + mapping of reads

- 5' clipping of amplification primers (optional)
- 3' removal of
  - Illumina adapters
  - low complexity stretches
  - low quality stretches
- filtering of non-SARS-CoV-2 reads (optional)
  - to our experience this step improves analyses quality a lot and hence we recommend that

 Temporary output:
 (will be removed at the end of the run)

- trimmed reads (fastq.gz)
- position sorted alignment file (bam)

## Step 2: variant calling + consensus sequence

- variant list (produced by Freebayes)
  - filtered for good coverage
  - filtered for >90% allele fraction
  - filtered against <10% allele fraction
- consensus sequences based on reference sequence
  - (examples given for variants in format Reference/Alternative:Depth\_Reference/Depth\_Alternative (a) G/C:0/80, (b) A/C:5/95, (c) C/T:50/50, (d) A/-:72/0, (e) G/-:10/0, (f) A/G:100/5 )
  - Variants detected in areas with low coverage (default: <20x) and poor quality (default: <10) are filtered out (i.e. N nucleotide reported)
  - (1) includes all variants with good coverage (default: >20x) and >90% allele fraction, all variants with <90% allele fraction are masked with N, regions with insufficient sequence depth are masked with N -- output suffix "masked_consensus.fasta"
    - examples (a) C, (b) C, (c) N, (d) A, (e) G, (f) A
  - (2) in addition to the one above including all mixed variants with <90% and >10% as IUPAC ambiguous bases -- output suffix: "iupac_consensus.fasta"
    - examples (a) C, (b) C, (c) Y, (d) A, (e) G, (f) A


## Step 2.1 variant annotation

If you want to get the variants annotated for their protein coding effect, please indicate so in the config file.
The tool used for annotation is SNPeff, which supports a variety of genomes.
To make sure that annotation works smoothly, please make sure that your reference genome is listed in the SNPeff supported genomes (Appendix).

## Step 3: generate summary files

- count mapped reads
- check reference base coverage
- investigate fragment size

Output:

- coverage
- frag_size.tsv (generated with an awk line) distribution of fragment sizes


## Step 4: generate report

This step may take a while in the first run of the pipeline as a bunch of packages have to be installed.

- optional runID is reported
- conditional table that warns the user, if samples that were identified as negative controls show high reference genome coverage
- table of read properties:
  - number of bases (before / after trimming)
    - if amplification primer clipping was done, be aware that "before trimming" does NOT show raw bases, but base count after primer clipping
  - length of reads (before / after trimming)
    - if amplification primer clipping was done, be aware that "before trimming" does NOT show raw read lengths, but read length after primer clipping
  - number of bases mapped (Q = 30)
- optional table listing the species filtering results emitted by Kraken
- table of mapping properties:
  - reads mapped to reference genome (number & fraction of input)
  - median / sd  of fragment size

- genome wide plot of coverage
- histogram of fragment sizes
- histogram / CDF of coverage vs number of bases (not implemented yet)

- table of reference coverage characteristics

Samples showing at least 20x sequencing depth at more than 95% of the reference genome are designated a successful genome sequencing.


## Troubleshoot

Please visit the [project's wiki ](https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe/-/wikis/home) for further information.

## Do you feel like contributing?

This [link](https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe/-/blob/edit_readme/CONTRIBUTING.md) will take you to our contribution guide

## Appendix

### SNPeff supported genomes

Listed are SARS-CoV-2 genomes supported by SNPeff.
Thanks to the Bioconda Team!

AP006557.1	SARS coronavirus TWH genomic RNA, complete genome.

AP006558.1	SARS coronavirus TWJ genomic RNA, complete genome.

AP006559.1	SARS coronavirus TWK genomic RNA, complete genome.

AP006560.1	SARS coronavirus TWS genomic RNA, complete genome.

AP006561.1	SARS coronavirus TWY genomic RNA, complete genome.

AY278488.2	SARS coronavirus BJ01, complete genome.

AY278489.2	SARS coronavirus GD01, complete genome.

AY278554.2	SARS coronavirus CUHK-W1, complete genome.

AY278741.1	SARS coronavirus Urbani, complete genome.

AY279354.2	SARS coronavirus BJ04, complete genome.

AY282752.2	SARS coronavirus CUHK-Su10, complete genome.

AY291451.1	SARS coronavirus TW1, complete genome.

AY304486.1	SARS coronavirus SZ3, complete genome.

AY304488.1	SARS coronavirus SZ16, complete genome.

AY304495.1	SARS coronavirus GZ50, complete genome.

AY310120.1	SARS coronavirus FRA, complete genome.

AY313906.1	SARS coronavirus GD69, complete genome.

AY323977.2	SARS coronavirus HSR 1, complete genome.

AY338174.1	SARS coronavirus Taiwan TC1, complete genome.

AY338175.1	SARS coronavirus Taiwan TC2, complete genome.

AY340092.1	SARS coronavirus BJ2232 RNA polymerase gene, partial cds.

AY345986.1	SARS coronavirus CUHK-AG01, complete genome.

AY345987.1	SARS coronavirus CUHK-AG02, complete genome.

AY345988.1	SARS coronavirus CUHK-AG03, complete genome.

AY348314.1	SARS coronavirus Taiwan TC3, complete genome.

AY350750.1	SARS coronavirus PUMC01, complete genome.

AY357075.1	SARS coronavirus PUMC02, complete genome.

AY357076.1	SARS coronavirus PUMC03, complete genome.

AY390556.1	SARS coronavirus GZ02, complete genome.

AY394978.1	SARS coronavirus GZ-B, complete genome.

AY394979.1	SARS coronavirus GZ-C, complete genome.

AY427439.1	SARS coronavirus AS, complete genome.

AY463059.1	SARS coronavirus ShanghaiQXC1, complete genome.

AY463060.1	SARS coronavirus ShanghaiQXC2, complete genome.

AY485277.1	SARS coronavirus Sino1-11, complete genome.

AY485278.1	SARS coronavirus Sino3-11, complete genome.

AY502923.1	SARS coronavirus TW10, complete genome.

AY502924.1	SARS coronavirus TW11, complete genome.

AY502925.1	SARS coronavirus TW2, complete genome.

AY502926.1	SARS coronavirus TW3, complete genome.

AY502927.1	SARS coronavirus TW4, complete genome.

AY502928.1	SARS coronavirus TW5, complete genome.

AY502929.1	SARS coronavirus TW6, complete genome.

AY502930.1	SARS coronavirus TW7, complete genome.

AY502931.1	SARS coronavirus TW8, complete genome.

AY502932.1	SARS coronavirus TW9, complete genome.

AY508724.1	SARS coronavirus NS-1, complete genome.

AY515512.1	SARS coronavirus HC/SZ/61/03, complete genome.

AY568539.1	SARS coronavirus GZ0401, complete genome.

AY572034.1	SARS coronavirus civet007, complete genome.

AY572035.1	SARS coronavirus civet010, complete genome.

AY572038.1	SARS coronavirus civet020, complete genome.

AY595412.1	SARS coronavirus LLJ-2004, complete genome.

AY613947.1	SARS coronavirus GZ0402, complete genome.

AY654624.1	SARS coronavirus TJF, complete genome.

AY686863.1	SARS coronavirus A022, complete genome.

AY686864.1	SARS coronavirus B039, complete genome.

AY687354.1	SARS coronavirus A001 spike glycoprotein gene, complete cds.

AY687355.1	SARS coronavirus A013 spike glycoprotein gene, complete cds.

AY687356.1	SARS coronavirus A021 spike glycoprotein gene, complete cds.

AY687357.1	SARS coronavirus A030 spike glycoprotein gene, complete cds.

AY687358.1	SARS coronavirus A031 spike glycoprotein gene, complete cds.

AY687359.1	SARS coronavirus B012 spike glycoprotein gene, complete cds.

AY687360.1	SARS coronavirus B024 spike glycoprotein gene, complete cds.

AY687361.1	SARS coronavirus B029 spike glycoprotein gene, complete cds.

AY687362.1	SARS coronavirus B033 spike glycoprotein gene, complete cds.

AY687364.1	SARS coronavirus B040 spike glycoprotein gene, complete cds.

AY714217.1	SARS Coronavirus CDC#200301157, complete genome.

AY772062.1	SARS coronavirus WH20, complete genome.

AY864805.1	SARS coronavirus BJ162, complete genome.

AY864806.1	SARS coronavirus BJ202, complete genome.

DQ022305.2	Bat SARS coronavirus HKU3-1, complete genome.

DQ084199.1	bat SARS coronavirus HKU3-2, complete genome.

DQ084200.1	bat SARS coronavirus HKU3-3, complete genome.

DQ182595.1	SARS coronavirus ZJ0301 from China, complete genome.

DQ640652.1	SARS coronavirus GDH-BJH01, complete genome.

DQ648856.1	Bat coronavirus (BtCoV/273/2005), complete genome.

DQ648857.1	Bat coronavirus (BtCoV/279/2005), complete genome.

EU371559.1	SARS coronavirus ZJ02, complete genome.

EU371560.1	SARS coronavirus BJ182a, complete genome.

EU371561.1	SARS coronavirus BJ182b, complete genome.

EU371562.1	SARS coronavirus BJ182-4, complete genome.

EU371563.1	SARS coronavirus BJ182-8, complete genome.

EU371564.1	SARS coronavirus BJ182-12, complete genome.

FJ882963.1	SARS coronavirus P2, complete genome.

GQ153539.1	Bat SARS coronavirus HKU3-4, complete genome.

GQ153540.1	Bat SARS coronavirus HKU3-5, complete genome.

GQ153541.1	Bat SARS coronavirus HKU3-6, complete genome.

GQ153542.1	Bat SARS coronavirus HKU3-7, complete genome.

GQ153543.1	Bat SARS coronavirus HKU3-8, complete genome.

GQ153544.1	Bat SARS coronavirus HKU3-9, complete genome.

GQ153545.1	Bat SARS coronavirus HKU3-10, complete genome.

GQ153546.1	Bat SARS coronavirus HKU3-11, complete genome.

GQ153547.1	Bat SARS coronavirus HKU3-12, complete genome.

JQ316196.1	SARS coronavirus HKU-39849 isolate UOB, complete genome.

JX993987.1	Bat coronavirus Rp/Shaanxi2011, complete genome.

JX993988.1	Bat coronavirus Cp/Yunnan2011, complete genome.

KC881005.1	Bat SARS-like coronavirus RsSHC014, complete genome.

KC881006.1	Bat SARS-like coronavirus Rs3367, complete genome.

KF367457.1	Bat SARS-like coronavirus WIV1, complete genome.

KF569996.1	Rhinolophus affinis coronavirus isolate LYRa11, complete genome.

KP886808.1	Bat SARS-like coronavirus YNLF_31C, complete genome.

KP886809.1	Bat SARS-like coronavirus YNLF_34C, complete genome.

MK062179.1	SARS coronavirus Urbani isolate icSARS, complete genome.

MK062180.1	SARS coronavirus Urbani isolate icSARS-MA, complete genome.

MK062181.1	SARS coronavirus Urbani isolate icSARS-C3, complete genome.

MK062182.1	SARS coronavirus Urbani isolate icSARS-C3-MA, complete genome.

MK062183.1	SARS coronavirus Urbani isolate icSARS-C7, complete genome.

MK062184.1	SARS coronavirus Urbani isolate icSARS-C7-MA, complete genome.

MN996532.1	Bat coronavirus RaTG13, complete genome.

NC_045512.2                                                 	COVID19 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1
