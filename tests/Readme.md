# Testing the Pipeline

[[_TOC_]]

This subfolder contains helper scripts and tooling, that allowe to validate the
pipeline. It is also used by the continous integration setup to run small tests
that verify the most basic functionality of the pipeline. 

## Quicktests

These are the small pipeline test runs that are mostly done to check syntax
and basic function of the pipeline

### [quicktest.sh](https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe/-/blob/edit_readme/tests/quicktest.sh)

| Test Requirement |  |
| ------ | ------ |
| bash | version 4 | 
| conda installation | 4.6.0+ |

From the project root simply run 

    >$ ./tests/quicktest.sh

## Comprehensive Quality Assurance



## Generation of Test Data

This section of the code base is a bit obscure.
It morphed from a small number of bash scripts into somthing akin to a 
workflow system. 

The basic pitch: 
    Do not store files. Only keep csv files with instructions how to generate 
    files you want. 

## Structure of the TestData Folder and how it works    

TL:DR:

    TestData/GEN_SAMPLES/gen_samples.sh Sample_7 gen

<blockquote>
Tries to create a 20/80 mixture of non modified and havily mutated Sars2 reads
 and highly contaminates the sample with a common cold corona virus in a ratio 11/89

Sample_2 is listed as the 20/80 mixture in the Gen_Samples folder and the same script
is called recursively on Sample_2.
Scripts calling scripts, leads to the Commoncold being downloaded from NCBI, 
aswell as the Wuhan seafood market Sars2, and the latter is also modified, 
reads are generated on all three sequences and mixed in the right ratios.  
</blockquote>


All get scipts work in the followig manner.

    get_script <regex> [info(default), cmd, run, generate]

**info** - Simply output the output files as absolute paths
This is the default execution mode

**cmd** - Output the commandlines it would run if it were in the run or generate mode

**run** - Run as usual and assume input files exist. 
All in script executions of other get_scipts are executed in info mode

**gen** -Force other getscripts to generate specified input files, if they do not exist. 

## Subfolders and Files 

### mutations.yaml

This file acts as a database of mutations that can be applyied 

### VCF_FILES / *get_vcfs.sh*
Generates custom vcf files based on mutations.yaml

### References
Downloads and aliases References from ncbi.
Uses Biopython and NCBIs Entrez-Utils to also allow search for Genbank and RefseQ
Assemblies based, when given an unspecific query.  
### MODIFIED_REFERENCE
Generates Modified reference fasta based on vcf and reference
### CONTAMINENTS
Simulates contamination reads (and unmodified reads)
### TARGET_RAW
Simulates target reads based on modified references
### GEN_SAMPLES
Creats the mixtures of reads and contaminants and illuminafies the simulated read
names










