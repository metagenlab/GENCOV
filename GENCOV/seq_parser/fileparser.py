#!/usr/bin/env python
"""fileparser script exemplifying how to use input utils

Executable script with main(cmd) as entry point
"""
import argparse
import os
import json
import sys
import re
import yaml
import input_utils


# Get path directory name of the script
try:
    script_path = os.path.dirname(os.path.realpath(__file__))
except:
    # Most likely life interpreter session __file__ not defined
    # Set to working directory
    input_utils.eprint("Does not properly run in interpreter session.\n"
                       'please strategically place a "import pdb; pdb.set_trace()"\n'
                        "Will set script_path to current working directory")
    script_path = os.getcwd()

def main(cmd=None):
    """ Entry Point for script

    Kwargs:
        CMD: default None (Allows to pass in cmd line)
    """
    parser = argparse_setup() # Setup Argument Parsers
    args = parser.parse_args(cmd) # Parse Arguments
    parse_files(args) # Run main Analysis


def argparse_setup():
    """Setup function for argument parser

    Encapsulates argument parser setup
    """
    parser = argparse.ArgumentParser() # Get Argument parser object
    parser.add_argument(
        '--input', '-i', dest='input',
        help=("input sample folder. Illumina filenames should end with"
              "_<lane>_<R1|R2>_number, e.g. Sample_12_345_R1_001.fastq,"
              " to find the right paired set."),
        required=True, nargs="+")
    parser.add_argument(
        '--sample_list', dest='sample_list',
        help=('tab separated sample list that allows entry of meta data.'
              'The file should be a tab separated list'
              ' with each line containing one sample.'
              ' The schema is: READ1 READ2 METADATA '
              ' Example: CHICG-134_S1_L001_R1_001<tab>'
              'CHICG-134_S1_L001_R2_001<tab>'
              'somenote1=interesting<tab>group=treatment'),
        required=False)
    parser.add_argument("-o", "--output", default=os.getcwd(), type=mkdir_if_not_exists)
    parser.add_argument("--outfile_name", default="samples.yaml")
    parser.add_argument("--fileguide",
                        default=os.path.join(script_path, 'filenames.yaml'))

    return parser # Return Parser



def parse_files(args):
    """Function doing all the work

    Parse Illumina and PacBio Sequences and Samples from file and
    directory list.

    Args:
        args: Namespace object returned by argument parser.
              Needs .input, .sample_list,
                    .fileguide, .output, .outfile_name

    """
    all_files = get_all_files(args.input)
    sample_list = args.sample_list
    filename_guide_file  = args.fileguide
    outfile_name = os.path.join(args.output, args.outfile_name)
    # Classify Files in input directories or sample list file
    formats_found,  discarded = input_utils.parse_sample_info_new(
            all_files, filename_guide_file, ['pacbio', 'illumina_fastq'],
            sample_list_file=sample_list)
    #Try to extract Illumina Data
    try:
        illumina_data = formats_found['il^lumina_fastq']
        # print(repr(format_known.mfrs).replace('>,', '>,\n'))
    except KeyError:
            exit('No samples found or none met criteria!!\n'
                 'These files were discarded:\n'
                 '%s' % '\n'.join(discarded))
    # try to extract pacbio data
    try:
        pacbio_data = formats_found['pacbio']
        print('looking for pacbio data...')
        pac_samples = pacbio_data.get_samples()
        with open(outfile_name.replace('.yaml',
                                           '_pacbio.yaml'),
                  'w') as sample_file:
            yaml.dump(pac_samples, sample_file, default_flow_style=False)
        print('looking for illumina data...')
    except KeyError:
        pass
    sample_dict = illumina_data.flatten_naive()
    # flatten naive just drops read info from sampel name
    # if it is paired end data.
    len_known = len(sample_dict)
    if False:
        try:
            # Last ditch effort to just get all the sequence files, that
            # don't match any standard and rename them.  Sticking them at the
            # end of the sample yaml
            salvaged_dict = illumina_data.leftovers.process_leftovers(
                rename=True,
                rename_start_index=len_known+1)
            sample_dict.update(salvaged_dict)
        except input_utils.AmbigiousPairedReadsError as err:
            input_utils.eprint('Failed parsing files with unrecognized'
                               ' naming convention\n',
                               'Reason:\n', err)
    # Write samples to working directory
    with open(outfile_name, 'w') as sample_file:
        yaml.dump(sample_dict, sample_file, default_flow_style=False)

def get_all_files(input_files):
    """ Get all files from a input file list

    Args:
        input_files:
    """
    results = []
    for file_or_dir in input_files:
        if os.path.isdir(file_or_dir):
            # Use walks in directories to traverse directory tree and
            # collect files
            for root, dirs, files in os.walk(file_or_dir):
                for _file in files:
                    if os.path.getsize(os.path.join(root, _file)) != 0:
                        results.append(os.path.join(root, _file))
        elif os.path.isfile(file_or_dir):
            if os.path.getsize(file_or_dir):
                results.append(file_or_dir)
    return results

def mkdir_if_not_exists(path):
    """Make directory if it does not exist

    As stated above, but also able to handle rare OSError race condition

    Args:
        path (:obj:`str`): path to directory to crate if it does not exist
    Raises:
        OSError: Most likely due to some other app creating the directory
        in between the isdir and mkdir call, leading to a race condition
        for the two processes.
    """
    if not os.path.isdir(path):
        try:
            os.makedirs(path)
            return path
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise
    return path

# When running as executable, auto call entry function main
if __name__ == "__main__":
    main()
