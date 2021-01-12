# -*- coding: utf-8 -*-
'''utils module for snakemake pipeline

This module provides utility classes and functions that make it easy to
handle nucliotide sequence data
'''

from itertools import takewhile
from collections import OrderedDict
import csv
import sys
import os.path
import operator
import re
import yaml
import pprint

''' FastQ Illumina filenames
FastQ with complex naming convention):
    <openBIS_Sequencing_Sample_Name>
    _<FlowCell_Name>_<Lane>_<Sample_Name>_<Index>_<sample number*>
    _<Lane>_<Read**>_<Running_Number***>_<Mismatches in Index****>
        openBIS ID contains underscores and "mismatches in index"
        field may or may not be present
Example:
    <BSSE_QGF_68928>_<CB94KANXX>_<6>_<Peteph_LPE5>_<TCTCGCGC_ATAGAGGC>
    _<S10>_<L006>_<R1>_<001>_<MM_1>.fastq.gz

    ([(openBIS_ID,3,'_'), flowcell_name, lane,
      (sample_name,['USER,ID'],'_'), sample_number]
    <BSSE>_<QGF>_<68928>_<CB94KANXX>_<6>_<Peteph_LPE5>_<TCTCGCGC_ATAGAGGC>
    _<S10>_<L006>_<R1>_<001>_<MM_1>.fastq.gz
(?P<openBIS_Sequencing_Sample_Name>)_(?P<FlowCell_Name>)_(?P<LaneSh>)_(?P<Sample_Name>)_(?P<Index>)_(?P<sample_number)*>_(?P<Lane>)_(?P<Read**>)_(?P<Running_Number>)_(?P<MismatchesIndex>).fastq.gz


IlluminaStandard:
    <Sample_Name>_<sample_number>_<Lane>_<Read**>_<Running_Number>
'''
# Exceptions

# illumina_standard_filename = [
#     ('sample_name', {}), ('sample_num', {'regex': 'S\d+'}),
#     ('lane', {'regex': 'L\d+'}), ('read', {'regex': 'R\d+'}),
#     ('running_num', {'regex': '\d+'})]


class SequenceFile:
    filename = ''
    path = ''
    ext = ''
    ext2 = ''
    add_exts = []

    def __init__(self, filename, path='', ext='', ext2='', add_exts=[]):
        self.filename = filename
        self.path = path
        self.ext = ext
        self.ext2 = ext2
        self.add_exts = []

    def __lt__(self, other):
        return (([self.filename, self.ext, self.ext2] + other)
                < ([other.filename, other.ext, other.ext2] + other.other))

    def __str__(self):
        return ''.join([os.path.join(self.path, self.filename),
                        self.ext, self.ext2] + self.add_exts)

    def __repr__(self):
        return ''.join([os.path.join(self.path, self.filename),
                        self.ext, self.ext2] + self.add_exts)


class Field():
    name = ''
    regex = ''
    subf_num = 1
    subf_sep = None
    optional = False

    def __init__(self, name, regex='', subf_num=1,
                 subf_sep=None, optional=False):
        self.name = name
        self.regex = regex
        self.subf_num = subf_num
        self.subf_sep = subf_sep
        self.optional = optional

    def to_regex(self, field_sep):
        if self.subf_num < 1:
            name = self.name
            raise MalformedFieldError(
                    'Field %s has less than 1 entry' % (name))
        elif self.subf_num == 1:
            field_regex = (r'[^%s]*' % field_sep if not self.regex
                           else self.regex)
            re_pattern = r'(?P<%s>%s)' % (self.name, field_regex)
        else:
            sfs = (field_sep
                   + (self.subf_sep if self.subf_sep is not None else ''))
            subf_sep = field_sep if self.subf_sep is None else self.subf_sep
            subf_regex = r'[^%s]*' % format(sfs)
            subf_regex = subf_sep.join([subf_regex]*self.subf_num)
            re_pattern = (
                    '(?P<%s>%s)' % (self.name, subf_regex))
        return re_pattern


class MalformedFieldError(Exception):
    '''Exception thrown if regex field is malformed
    '''
    pass


class AmbigiousPairedReadsError(Exception):
    ''' Exception thrown if read pairs can't be matched

    Paired Read data has at most two files that are grouped together
    If more are found, this esxception is raised.
    Example:
        S1_R1.fastq, S1_R2.fastq S1_R3.fastq > Exception
        S1_R1.fastq, S1_R2.fastq S1_R2.fastq > Exception
        S1_R1.fastq, S1_R2.fastq  > OK
    '''
    pass


class UnknownExtensionError(Exception):
    '''Exception thrown if input sequence file has an unknown file extension

    If the pipeline only excepts nucleotide sequence files, like
    Fastas and Fastqs, inputting a mapping file or index will
    cause this exception to be thrown.
    '''
    pass


class IncongruentFieldValueError(Exception):
    '''Exception thrown if ther is a missmatch in a sample file grouping

    If a grouping of files, that belong to the same SampleId grouping, differ
    in a non variable field of the filename, this exception is raised
    '''
    pass


class SampleIDNotUniqueError(Exception):
    '''File exists more than once in sample file grouping

    This might be caused by a sample id that is not unique
    '''
    pass


class FormatMismatch(Exception):
    pass


class SanityCheckFailedError(Exception):
    pass

# Classes


class SampleManager():
    curr_sample_dict = {}
    main_sample_dict = {}
    selected_sample_collection = None
    sample_collection_labels = set()
    parsed_technologies = set()
    enforce_type_consistency = True
    parsed_MFR_Collections = {}

    def __check_sample_dict_valid(self, sample_dict):
        return all(sample_dict.keys()) and all(
                issubclass(type(sample_dict[x]), SampleInfo)
                for x in sample_dict)

    def __check_not_unique_and_reason(self, dict1, dict2):
        output = []
        for sample in dict1:
            if sample in dict2:
                output.append(
                        'Collision:\n%s\n\t --- \t\n%s\n' % (
                            repr(dict1[sample]),
                            repr(dict2[sample])))
        return output

    def __init__(self):
        self.curr_sample_dict = dict()
        self.main_sample_dict = dict()
        self.selected_sample_collection = 'default'
        self.sample_collection_labels = set()
        self.parsed_technologies = set()
        self.enforce_type_consistency = True

    def has_samples(self):
        return bool(self.curr_sample_dict or self.main_sample_dict)

    def has_samples_selected(self):
        return bool(self.curr_sample_dict)

    def load_from_sample_dict(self, sample_dict, overwrite=False):
        valid = self.__check_sample_dict_valid(sample_dict)
        if not valid:
            raise Exception("Sample dict not valid:"
                            " Cannot have empty ids or values that do not\n"
                            " belong to a subclass of input_utils.SampleInfo")
        if not overwrite:
            not_unique_and_reason = self.__check_not_unique_and_reason(
                    sample_dict, self.curr_sample_dict)
            if not_unique_and_reason:
                raise Exception(
                        'Error: SampleID collision\n while combining'
                        ' loaded sample dict and current sample dict\n'
                        'err:\n %s' % '\n'.join(not_unique_and_reason))
            self.curr_sample_dict.update(sample_dict)

    def load_from_sample_dict_file(self, sample_dict_file, overwrite=False):
        is_path = isinstance(sample_dict_file, str)
        if is_path:
            with open(sample_dict_file) as s_fh:
                sample_dict = yaml.load(s_fh)
        else:
            sample_dict = yaml.load(sample_dict_file)
        valid = self.__check_sample_dict_valid(sample_dict)
        if not valid:
            raise Exception("Sample dict not valid:"
                            " Cannot have empty ids or values that do not\n"
                            " belong to a subclass of input_utils.SampleInfo")
        if not overwrite:
            not_unique_and_reason = self. __check_not_unique_and_reason(
                    sample_dict, self.curr_sample_dict)
            if not_unique_and_reason:
                raise Exception(
                        'Error: SampleID collision\n while combining'
                        ' loaded sample dict and current sample dict\n'
                        'err:\n %s' % '\n'.join(not_unique_and_reason))
            self.curr_sample_dict.update(sample_dict)

    def sample_dict_to_tsv(self, out_file='samplelist.tsv', target_label=''):
        meta_header = self.get_meta_field_names(target_label)
        if not target_label:
            target = self.curr_sample_dict
            target_label = self.selected_sample_collection
        else:
            target = self.main_sample_dict[target_label]

        if self.target_label in ['default', 'illumina']:
            Sequence_files = OrderedDict([
                    ('READ1', (lambda x: x.READ1)),
                    ('READ2', (lambda x: ''
                               if not hasattr(x, 'READ2')
                               else x.READ2))])
        sequence_group = list(Sequence_files.keys())
        with open(out_file, 'w') as out_fh:
            tsv_writer = csv.writer(out_fh, delimiter='\t')
            tsv_writer.writerow(['ID'] + sequence_group + meta_header)
            for sample_id in target:
                row = ([sample_id]
                       + [Sequence_files[x](target[sample_id])
                          for x in Sequence_files])
                for info in meta_header:
                    try:
                        row.append(target[sample_id].meta_info[info])
                    except KeyError:
                        row.append('')
                tsv_writer.writerow(row)

    def parse_samples_from_tsv(self, file_input,
                               use_common_fields_as_id=False):
        pass

    def parse_samples_from_file_list(self, files,
                                     use_common_fields_as_id=False):
        pass

    def process_parsed_files(self, mode='', unknown_formats=False):
        pass

    def filter_samples(self, target_label='',
                       select_add_info=[], select_meta_info=[]):
        """Todo: think of new label for filtered
        """
        no_rename = False
        if not target_label:
            target = self.curr_sample_dict
            target_label = self.selected_sample_collection
            no_rename = True
        else:
            try:
                target = self.main_sample_dict[target_label]
            except KeyError as err:
                raise KeyError('%s not found in main_sample_dict\n'
                               'Error: %s' % (target_label, err))
        selected = set(x[0] for x in sample_dict_get_info(
                target, select_add_info=select_add_info,
                select_meta_info=select_meta_info))
        self.curr_sample_dict = dict((x, target[x]) for x in target
                                     if x in selected)
        if not no_rename:
            self.selected_sample_collection = target_label + '_filtered'

    def get_meta_field_names(self, target_label=''):
        field_names = set()
        if not target_label:
            target = self.curr_sample_dict
        else:
            try:
                target = self.main_sample_dict[target_label]
            except KeyError as err:
                raise KeyError('%s not found in main_sample_dict\n'
                               'Error: %s' % (target_label, err))
        for id_string in target:
            for meta_field_name in target[id_string].meta_info:
                field_names.add(meta_field_name)
        return sorted(list(field_names))


def sample_dict_get_info(sample_dict,
                         select_add_info=[], select_meta_info=[],
                         print_add_info=[], print_meta_info=[],
                         ommit_id=False):
    '''
    '''
    def __selector_valid(selector):
        operators = ['!', '==', '<=', '>=', '<', '>', '!=']
        formats = ['str', 'num', 'date']
        if isinstance(selector, tuple):
            try:
                f_num = len(selector)
                for x in selector:
                    __selector_valid(x)
                if selector[1] not in operators:
                    raise ValueError('Operator should be one of these\n'
                                     '%s' % ' '.join(operators))
                if f_num == 2 and selector[1] is not '!':
                    raise ValueError('selector tuple of size 2:\n'
                                     'operator needs to be not ("!")')
                if f_num < 2 and f_num > 4:
                    raise ValueError('selector tuple should'
                                     ' have 3 or 4 fields')
                if len(selector) == 4 and selector[3] not in formats:
                        raise ValueError('Expectet format strings are:\n'
                                         '%s' % ' '.join(formats))
            except ValueError as err:
                raise ValueError('TupleSelectorError:'
                                 ' It violated the format\n'
                                 '(<selector_str>, <operator_str>,'
                                 ' <value_str>,\n'
                                 ' optional |<str,num,date>|)\nerr:%s' % err)
        elif not isinstance(selector, str):
            raise ValueError('selectors')
        return selector

    def __print_selector_valid(selector):
        if not isinstance(selector, str):
            raise ValueError('Error: Print selector should be a string')
        return selector

    def __is_selection(sample, select_add_info, select_info):

        def __comparison(sample, selector,
                         access_fun=lambda sample, x: sample[x]):
            op_fun = {
                '<': operator.lt,
                '<=': operator.le,
                '==': operator.eq,
                '!=': operator.ne,
                '>=': operator.ge,
                '>': operator.gt,
                '!': operator.not_,
                'not': operator.not_}
            if isinstance(selector, tuple):
                selector = __selector_valid(selector)
                if len(selector) == 2:
                    return not access_fun(sample, selector[0])
                elif len(selector) == 4:
                    Exception('Format mode not implemented yet')
                else:
                    try:
                        return op_fun[selector[1]](*[
                            access_fun(sample, selector[0]), selector[2]])
                    except KeyError:
                        return False
            else:
                try:
                    access_fun(sample, selector)
                except KeyError:
                    return False
        # Check add_info
        if not all(__comparison(
                        sample, selector,
                        access_fun=lambda s, x: s.add_info[x])
                   for selector in select_add_info):
            return False
        # Check meta_info
        if not all(__comparison(
                        sample, selector,
                        access_fun=lambda s, x: s.meta_info[x])
                   for selector in select_meta_info):
            return False
        return True

    def __get_print_fields(sample, print_add_info, print_meta_info, ommit_id):
        out = [] if ommit_id else [sample.ID]
        out.extend((sample.add_info[__print_selector_valid(x)]
                    for x in print_add_info))
        out.extend((sample.meta_info[__print_selector_valid(x)]
                    for x in print_meta_info))
        return out

    return [__get_print_fields(sample_dict[s],
                               print_add_info,
                               print_meta_info, ommit_id)
            for s in sample_dict
            if __is_selection(sample_dict[s],
                              select_add_info, select_meta_info)]


class MultiFileReads:
    def __init__(self):
        pass

    def get_files(self):
        pass

    def add_files(self):
        pass


class IlluminaMFR(MultiFileReads):
    n_files = 0
    uses_common_fields_as_id = False
    n_branches = 0
    format_name = ''
    field_names = []
    id_fields = []
    id_string = ''
    id_sep = ''
    immutable_id = False
    main_sep = ''
    meta_info = {}
    var_fields = []
    ignore_fields = []
    regex_str = ''
    regex = None
    sample_dict = {}
    non_var_dict = {}
    common_fields_dict = {}
    default_field_values = {}

    def __init__(self, filename_rule, main_sep='_',
                 id_sep='', id_fields=['sample_name'],
                 immutable_id=False,
                 var_fields=['lane', 'running_num'],
                 format_name='new_format',
                 default_field_values={'sample_num': 'S1',
                                       'lane': 'L001',
                                       'read': 'R1',
                                       'running_num': '001'}):
        super().__init__()
        # Construct regex string
        self.regex_str = main_sep.join([
                Field(name, **key_opts).to_regex(main_sep)
                for (name, key_opts) in filename_rule])
        self.regex = re.compile(self.regex_str)
        self.field_names = [name for (name, _) in filename_rule]
        self.id_fields = id_fields
        self.id_sep = id_sep
        self.id_string = ''
        self.immutable_id = immutable_id
        self.default_field_values = default_field_values
        self.var_fields = var_fields
        self.sample_dict = {}
        self.non_var_dict = {}
        self.format_name = format_name
        self.main_sep = main_sep
        self.meta_dict = dict()
        self.n_files = 0
        self.n_branches = 0
        self.common_fields_dict = {}

    def __repr__(self):
        return 'IlluminaMFR<%s>' % ', '.join(
                ['format_name: %s' % self.format_name,
                 '\nn_files: %i' % self.n_files,
                 'n_branches: %i' % self.n_branches,
                 'tree:\n%s' % pprint.pformat(self.sample_dict)])

    def add_files(self, file_list, read_target=0):
        if isinstance(file_list, SequenceFile):
            file_list = [file_list]
        file_list = sorted(file_list)
        for file in file_list:
            try:
                field_dict = self.regex.match(file.filename).groupdict()
            except AttributeError:
                raise FormatMismatch(
                    'FormatMismatch: %s does not fit %s' % (file,
                                                            self.format_name))
            sample_dict = self.sample_dict
            read = None
            var_dict = {}
            if self.n_files == 0:
                self.common_fields_dict = dict((key, val) for (key, val)
                                               in field_dict.items()
                                               if key
                                               not in (self.var_fields
                                                       + ['read']))
                self.update_id_field()
            for field in self.field_names:
                f_value = field_dict[field]
                sample_key = field + ':' + f_value
                read_field = field == 'read'
                if read_field:
                    read_id = f_value[-1]
                    if not read_target or read_id == str(read_target):
                        read = file.filename
                elif field in self.var_fields:
                    if self.n_files == 0 or sample_key not in sample_dict:
                        sample_dict[sample_key] = {}
                    sample_dict = sample_dict[sample_key]
                    var_dict[field] = f_value
                    # print(var_dict)
                else:
                    if self.n_files == 0:
                        self.non_var_dict[sample_key] = 42
                    elif (not read_field
                          and sample_key not in self.non_var_dict):
                        raise IncongruentFieldValueError(
                            'File: %s\n of ID Group %s has'
                            ' group missmatch in field: %s\n'
                            ' with value %s' % (
                                file, self.id_string, field, f_value))
            if read is None:
                read_id = 1
                read = file.filename
            # if read is not None: #
            if 'read%s' % read_id in sample_dict:
                raise SampleIDNotUniqueError(
                        'Error: files \n%s\n%s\n'
                        'have the same variable arguments. ' % (
                            str(file), sample_dict['read'+str(read_id)]))
            var_dict['read%s' % read_id] = file
            if len(sample_dict) < 1:
                self.n_branches += 1
            sample_dict.update(var_dict.copy())
            self.n_files += 1

    def update_id_field(self, id_sep=None, id_fields=None):
        if self.immutable_id:
            return
        if id_sep is not None:
            self.id_sep = id_sep
        if id_fields:
            self.id_fields = id_fields
        self.id_string = self.id_sep.join(
                [self.common_fields_dict[field] for field in self.id_fields])
        if id_sep or id_fields:
            self.use_common_fields_as_id = False

    def set_format_name(self, new_name):
        self.format_name = new_name

    def get_files(self, undetermined=None):
        output = []
        to_do = [self.sample_dict]
        while to_do:
            curr_node = to_do.pop()
            if not any(isinstance(next_node, dict)
                       for next_node in curr_node.values()):
                output += [curr_node]
            else:
                to_do.extend([curr_node[n] for n in curr_node])
        return output

    def set_meta_info(self, meta_info_dict):
        self.meta_info = meta_info_dict

    def set_immutable_id(self, id, add_standard_suffix=False):
        self.immutable_id = True
        self.id_string = id + ('S1_L001_001' if add_standard_suffix else '')

    def use_common_fields_as_id(self):
        new_id_fields = [f for f in self.field_names
                         if f not in (self.var_fields+['read'])]
        self.update_id_field(id_fields=new_id_fields)
        self.uses_common_fields_as_id = True


class PacBioMFR(MultiFileReads):
    n_files = 0
    n_branches = 0
    format_name = ''
    field_names = []
    id_fields = []
    id_string = ''
    id_sep = ''
    main_sep = ''
    var_fields = []
    ignore_fields = []
    regex_str = ''
    regex = None
    sample_dict = {}
    non_var_dict = {}
    common_fields_dict = {}

    def __init__(self, filename_rule, main_sep='_',
                 id_sep='', id_fields=['sample_name'],
                 var_fields=['dummy'],
                 format_name='new_format'):
        super().__init__()
        # Construct regex string
        self.regex_str = main_sep.join([
                Field(name, **key_opts).to_regex(main_sep)
                for (name, key_opts) in filename_rule])
        self.regex = re.compile(self.regex_str)
        self.field_names = [name for (name, _) in filename_rule]
        self.id_fields = id_fields
        self.id_sep = id_sep
        self.id_string = ''
        self.var_fields = var_fields
        self.sample_dict = {}
        self.non_var_dict = {}
        self.format_name = format_name
        self.main_sep = main_sep
        self.n_files = 0
        self.n_branches = 0
        self.common_fields_dict = {}

    def __repr__(self):
        return 'PacBioMFR<%s>' % ', '.join(
                ['format_name: %s' % self.format_name,
                 '\nn_files: %i' % self.n_files,
                 'n_branches: %i' % self.n_branches,
                 'tree:\n%s' % pprint.pformat(self.sample_dict)])

    def add_files(self, file_list, file_selector_target=''):
        if isinstance(file_list, SequenceFile):
            file_list = [file_list]
        file_list = sorted(file_list)
        for file in file_list:
            try:
                field_dict = self.regex.match(file.filename).groupdict()
            except AttributeError:
                raise FormatMismatch(
                    'FormatMismatch: %s does not fit %s' % (file,
                                                            self.format_name))
            sample_dict = self.sample_dict
            file_selector = None
            var_dict = {}
            if self.n_files == 0:
                self.common_fields_dict = dict((key, val) for (key, val)
                                               in field_dict.items()
                                               if key
                                               not in (self.var_fields
                                                       + ['file_selector']))
                self.update_id_field()
            for field in self.field_names:
                f_value = field_dict[field]
                sample_key = field + ':' + f_value
                file_selector_field = field == 'file_selector'
                if file_selector_field:
                    file_selector_id = f_value
                    # Detect if we are in the variable filename field
                    in_field_that_determines_file_identity = (
                            not file_selector_target
                            or file_selector_id == str(file_selector_target))
                    if in_field_that_determines_file_identity:
                            file_selector = file.filename
                elif field in self.var_fields:
                    if self.n_files == 0 or sample_key not in sample_dict:
                        sample_dict[sample_key] = {}
                    sample_dict = sample_dict[sample_key]
                    var_dict[field] = f_value
                else:
                    if self.n_files == 0:
                        self.non_var_dict[sample_key] = 42
                    elif (not file_selector_field
                          and sample_key not in self.non_var_dict):
                        raise IncongruentFieldValueError(
                            'File: %s\n of ID Group %s has'
                            ' group missmatch in field: %s\n'
                            ' with value %s' % (
                                file, self.id_string, field, f_value))
            if 'file_selector_id' in sample_dict:
                raise SampleIDNotUniqueError(
                        'Error: files \n%s\n%s\n'
                        'have the same variable arguments. ' % (
                            file_selector,
                            sample_dict[str(file_selector_id)]))
            if file_selector is not None:
                var_dict[file_selector_id] = file
            if len(sample_dict) < 1:
                self.n_branches += 1
            sample_dict.update(var_dict.copy())
            self.n_files += 1

    def update_id_field(self, id_sep=None, id_fields=None):
        if id_sep is not None:
            self.id_sep = id_sep
        if id_fields:
            self.id_fields = id_fields
        self.id_string = self.id_sep.join(
                [self.common_fields_dict[field] for field in self.id_fields])

    def set_format_name(self, new_name):
        self.format_name = new_name

    def get_files(self, undetermined=None):
        output = []
        to_do = [self.sample_dict]
        while to_do:
            curr_node = to_do.pop()
            if not any(isinstance(next_node, dict)
                       for next_node in curr_node.values()):
                output += [curr_node]
            else:
                to_do.extend([curr_node[n] for n in curr_node])
        return output


class IonTorrentMFR(MultiFileReads):
    pass


def test_fun(filename):
    '''testfunction used for debugging purposes

    imports first format of given filenames.yaml, builds and returns
    corresponding multi file read object
    '''
    formats = get_formats_from_file(filename)
    basel = list(formats[0].values())[0]
    mfr = IlluminaMFR(basel['format'].items(),
                      main_sep=basel['main_sep'])
    return mfr


class SampleInfo:
    ID = ''
    add_info = {}
    meta_info = {}

    def __init__(self):
        pass


class SampleInfoPaired(SampleInfo):
    ''' Information Container for Paired End Read Data

    Attributes:
        ID (str): first section of file base name shared by read pair
        READ1 (str): Absolute path of first half of read pairs
        READ2 (str): Absolute path to second half of read pairs
        exts (:obj:`list` of :obj:`str`):
            List of extensions used by read pair in correct order
        zip_exts (:obj:`list` of :obj:`str`):
            List of compression extensions used by pair (also in correct order)
            If they are empty strings, no compression extension is used
    '''
    ID = ''
    READ1 = ''
    READ2 = ''
    read_ids = ['', '']
    exts = ['', '']
    zip_exts = ['', '']
    add_info = {}
    meta_info = {}

    def __init__(self, r1, r2, id_str, read_ids, exts=['.fastq', '.fastq'],
                 zip_exts=['', ''], add_info={}, meta_info={}):
        '''Initializer

        Reads in the attributes in the order
            read1, read2, ids, exts, zip, exts

        Kwargs:
            exts: default ['.fastq', '.fastq']
            zip_exts: default ['', '']
        '''
        self.ID = id_str
        self.READ1 = r1
        self.READ2 = r2
        self.exts = exts
        self.read_ids = read_ids
        self.zip_exts = zip_exts
        self.add_info = add_info
        self.meta_info = meta_info

    def __str__(self):
        return 'SampleInfoPaired()'

    def __repr__(self):
        att = (self.ID, self.READ1, self.READ2, ','.join(self.exts),
               ','.join(self.zip_exts))
        return ('SampleInfoPaired()'
                if not any(att)
                else '<SampleInfoPaired:'
                     ' ID=%s,\n READ1=%s,\n READ2=%s,\n'
                     ' exts=%s, zip_exts=%s>' % att)


class SampleInfoSingle(SampleInfo):
    ''' Information Container for Single End Read Data

    Attributes:
        ID (str): first section of read file base name
        READ1 (str): Absolute path to read file
        ext ( :obj:`str`):
            Extension used by read data file
        zip_ext (:obj:`str`):
            Compression extension used by read file
            If the string is empty, no compression extension is used
    '''
    ID = ''
    READ1 = ''
    exts = ''
    read_ids = ''
    zip_exts = ''
    add_info = {}
    meta_info = {}

    def __init__(self, ids, r1, read_ids,
                 ext='.fastq', zip_ext='', add_info={}, meta_info={}):
        '''Initializer

        Reads in the attributes in the order
            ids, read1, exts, zip, exts

        Kwargs:
            ext: default  '.fastq'
            zip_ext: default ''
        '''
        self.ID = ids
        self.READ1 = r1
        self.exts = ext
        self.read_ids = read_ids
        self.zip_exts = zip_ext
        self.add_info = add_info
        self.meta_info = meta_info

    def __str__(self):
        return 'SampleInfoSingle()'

    def __repr__(self):
        att = (self.ID, self.READ1, self.ext, self.zip_ext)
        return ('SampleInfoSingle()'
                if not any(att)
                else '<SampleInfoSingle:'
                     ' ID=%s,\n READ1=%s,\n ext=%s, zip_ext=%s>' % att)


class PacBioSampleInfoRS_II:
    bax1 = ''
    bax2 = ''
    bax3 = ''
    metadata = ''
    bas = ''
    add_info = {}

    def __init__(self, ids, metadata='', bax1='', bax2='', bax3='',
                 bas='', add_info={}):
        '''Initializer

        Reads in the attributes in the order
            ids, read1, exts, zip, exts

        Kwargs:
            ext: default  '.fastq'
            zip_ext: default ''
        '''
        self.ID = ids
        self.metadata = metadata
        self.bas = bas
        self.add_info = add_info
        self.bax1 = bax1
        self.bax2 = bax2
        self.bax3 = bax3

    def __str__(self):
        return 'SampleInfoSingle()'


class PacBioSampleInfoRS:
    ID = ''
    metadata = ''
    bas = ''
    add_info = {}

    def __init__(self, ids, metadata='', bas='', add_info={}):
        '''Initializer

        Reads in the attributes in the order
            ids, read1, exts, zip, exts

        Kwargs:
            ext: default  '.fastq'
            zip_ext: default ''
        '''
        self.ID = ids
        self.metadata = metadata
        self.bas = bas
        self.add_info = add_info

    def __str__(self):
        return 'SampleInfoSingle()'


class ReferenceInfo:
    ''' Information Container for Genomic Reference Data

    Attributes:
        ID (str): reference file base name without extension
        REFERENCE (str): Absolute path to reference file
        ext ( :obj:`str`):
            Extension used by reference data file
        zip_ext (:obj:`str`):
            Compression extension used by reference file
            If the string is empty, no compression extension is used
    '''
    ID = ''
    REFERENCE = ''
    ext = ''
    zip_ext = ''

    def __init__(self, id, reference, ext='.fna', zip_ext='.gz'):
        self.ID = id
        self.REFERENCE = reference
        self.ext = ext
        self.zip_ext = zip_ext

    def __str__(self):
        return 'ReferenceInfo()'

    def __repr__(self):
        att = (self.ID, self.REFERENCE, self.ext, self.zip_ext)
        return ('ReferenceInfo()'
                if not any(att)
                else '<ReferenceInfo:'
                     ' ID=%s,\n Reference=%s,\n ext=%s, zip_ext=%s>' % att)


def eprint(*args, **kwargs):
    '''
    print function that prints to stderr

    :return: returns nothing
    '''
    print(*args, file=sys.stderr, **kwargs)


def test_extension(filename, extension_list):
    ''' tests which extension a file uses

    Args:
        filename (:obj:`str`):
            name of file whose extension will get checked
        extension_list (:obj:`list` of :obj:`str`):
            list of extensions that the file will be checked against.
            should contain the dot and extensions that share a prefix should
            be sorted in ascending order

    Returns:
        (:obj:`str`): Extension used or '' if not found in extension_list
    '''
    res = ''
    for ext in extension_list:
        if len(filename.split(ext)) == 2:
            res = ext
            break
    return res


def parse_sample_info(sample_list, format_dict,
                      use_common_fields_as_id=False,
                      target_formats=['illumina_fastq']):
    '''Parses filenames and generates SampleInfoObjects

    Turns list of input files into read data containers.
    It finds pairs for paired end read data and determines
    sample ids, which compression extension is used (if any),
    and which nucleotide sequence file extension is used.

    It only accepts files that end in <seqf_ext>(otional:<comp_ext>)

    Args:
        sample_list(:obj:`list` of :obj:`str`): list of filenames
        format_dict(:obj:`dict` nested and with various types):
            file naming conventions loaded from yaml file

    Kwargs:
        use_common_fields_as_id(:obj:`bool` default False):
            Ignores ID Fields specified in formats config and just uses
            all fields that are not marked as variable as ID
        target_formats(:obj:`list` of :obj:`str`):
            List of machine target_formats to check against
            Choose from:
                illumina_fastq, ion_torrent_bam, pacbio
    Returns:
        :obj:`dict` of :obj:`MFR_Collection`
            Returns dictionary of MFR_Collections found,
            with respective target_format as key value.
            Example:
                result["illumina_fastq"] returns a Illumina_MFR_Collection()

            It also provides a dictionary entry for discarded files
                result["discarded"]

    Raises:
        UnknownExtensionError: If sequence file extension unknown
        AmbigiousPairedReadsError: If paired data has to many matching files
    '''
    collections = {'ion_torrent_bam': IonTorrent_MFR_Collection,
                   'illumina_fastq': Illumina_MFR_Collection,
                   'pacbio': PacBio_MFR_Collection}
    sample_list = sorted(sample_list)
    mfr_samples = {}  # Collection of detected multifile samples
    mfrs_found = dict()
    discarded = list()
    # Accumulate sample information and build samples dictinary
    for sample in sample_list:
        found_format_or_is_leftover = False
        for target_format_type in target_formats:
            if found_format_or_is_leftover:
                break
            seq_exts = format_dict[target_format_type]['main_exts']
            format_list = format_dict[target_format_type]['formats']
            # Select which mfr type to try
            mfr_type = collections[target_format_type].mfr_type
            try:
                zip_exts = format_dict[target_format_type]['secondary_exts']
            except KeyError:
                zip_exts = []

            used_ext = test_extension(sample, seq_exts)
            if not used_ext:
                continue
                raise UnknownExtensionError(
                        'Extension not recognized\n%s' % sample)
            sample_string, zipped = sample.split(used_ext)
            if zipped and not test_extension(zipped, zip_exts):
                continue
            path = os.path.dirname(sample)
            sample = os.path.basename(sample_string)
            # Get first section of file.filename for ID
            seq_file = SequenceFile(sample, path=path,
                                    ext=used_ext, ext2=zipped, add_exts=[])
            # Check if known format
            mfr = find_format(seq_file, format_list, mfr_type)
            if target_format_type not in mfrs_found:
                mfrs_found[target_format_type] = (
                        collections[target_format_type]())
            if issubclass(type(mfr), MultiFileReads):
                if use_common_fields_as_id:
                    mfr.use_common_fields_as_id()
                mfrs_found[target_format_type].add(mfr, seq_file)
            else:
                mfrs_found[target_format_type].leftovers.add(sample, path,
                                                             zipped, used_ext)
            found_format_or_is_leftover = True
        if not found_format_or_is_leftover:
            discarded.append(sample)

    # Process samples and build
    for id_str in mfr_samples:
        mfr = mfr_samples[id_str]
        print('MultiFileSample: ', id_str, 'format:', mfr.format_name)
        print(mfr.get_files())

    return mfrs_found, discarded


def parse_sample_info_new(input_file_list, format_guide_file,
                          use_common_fields_as_id=False,
                          target_formats=['illumina_fastq'],
                          sample_list_file=''):
    '''Parses filenames and generates SampleInfoObjects

    Turns list of input files into read data containers.
    It finds pairs for paired end read data and determines
    sample ids, which compression extension is used (if any),
    and which nucleotide sequence file extension is used.

    It only accepts files that end in <seqf_ext>(otional:<comp_ext>)

    Args:
        input_file_list(:obj:`list` of :obj:`str`): list of filenames
        format_guide_file(:obj:`str` path to yaml file with
                          file naming conventions):
            file naming conventions are loaded from this yaml file

    Kwargs:
        sample_list_file(:obj:`str` default: empty string):
            path/to/<sample_list_file>.tsv
        use_common_fields_as_id(:obj:`bool` default False):
            Ignores ID Fields specified in formats config and just uses
            all fields that are not marked as variable as ID
        target_formats(:obj:`list` of :obj:`str`):
            List of machine target_formats to check against
            Choose from:
                illumina_fastq, ion_torrent_bam, pacbio
    Returns:
        :obj:`dict` of :obj:`MFR_Collection`
            Returns dictionary of MFR_Collections found,
            with respective target_format as key value.
            Example:
                result["illumina_fastq"] returns a Illumina_MFR_Collection()

            It also provides a dictionary entry for discarded files
                result["discarded"]

    Raises:
        UnknownExtensionError: If sequence file extension unknown
        AmbigiousPairedReadsError: If paired data has to many matching files
    '''
    collections = {'ion_torrent_bam': IonTorrent_MFR_Collection,
                   'illumina_fastq': Illumina_MFR_Collection,
                   'pacbio': PacBio_MFR_Collection}
    # trying to read format guide file
    failed_reading_formats_file = ''
    try:
        format_dict = get_formats_from_file(format_guide_file)
    except yaml.YAMLError as err:
        failed_reading_formats_file = 'Error parsing yaml: %s' % err
    except IOError as err:
        failed_reading_formats_file = 'Error opening file: %s' % err
    except Exception as err:
        failed_reading_formats_file = 'Unknown_exception: %s' % err
    # Checking if tsv sample list with meta data was provided
    if failed_reading_formats_file:
        eprint('Warning: error encountered while reading filename guide file\n'
               '    %s' % failed_reading_formats_file)
    mfrs_found_tsv = dict()
    if sample_list_file:
        mfrs_found_tsv = _determine_mfr_type_from_sample_list(
                sample_list_file, collections, format_dict, target_formats,
                use_common_fields_as_id)
    # parse sample lists
    mfrs_found = dict()
    # Initialize format types that were read in the tsv sample sheet logic
    for format_type in mfrs_found_tsv:
        mfrs_found[format_type] = (
            collections[format_type]())
    input_file_list = sorted(input_file_list)
    mfr_samples = {}  # Collection of detected multifile samples
    discarded = list()
    # Accumulate sample information and build samples dictinary
    for sample in input_file_list:
        found_format_or_is_leftover = False
        for target_format_type in target_formats:
            if found_format_or_is_leftover:
                break
            seq_exts = format_dict[target_format_type]['main_exts']
            format_list = format_dict[target_format_type]['formats']
            # Select which mfr type to try
            mfr_type = collections[target_format_type].mfr_type
            try:
                zip_exts = format_dict[target_format_type]['secondary_exts']
            except KeyError:
                zip_exts = []

            used_ext = test_extension(sample, seq_exts)
            if not used_ext:
                continue
                raise UnknownExtensionError(
                        'Extension not recognized\n%s' % sample)
            sample_string, zipped = sample.split(used_ext)
            if zipped and not test_extension(zipped, zip_exts):
                continue
            path = os.path.dirname(sample)
            sample = os.path.basename(sample_string)
            # Get first section of file.filename for ID
            seq_file = SequenceFile(sample, path=path,
                                    ext=used_ext, ext2=zipped, add_exts=[])
            # Check if known format
            mfr = find_format(seq_file, format_list, mfr_type)
            if target_format_type not in mfrs_found:
                mfrs_found[target_format_type] = (
                        collections[target_format_type]())
                if target_format_type not in mfrs_found_tsv:
                    mfrs_found_tsv[target_format_type] = (
                            collections[target_format_type]())

            if issubclass(type(mfr), MultiFileReads):
                if use_common_fields_as_id:
                    mfr.use_common_fields_as_id()
                if mfr.id_string in mfrs_found_tsv[target_format_type].mfrs:
                    found_format_or_is_leftover = True
                    continue
                mfrs_found[target_format_type].add(mfr, seq_file)
            else:
                # ToDo: detect if leftover was already processed
                #       in sample list file logic
                mfrs_found[target_format_type].leftovers.add(sample, path,
                                                             zipped, used_ext)
            found_format_or_is_leftover = True
        if not found_format_or_is_leftover:
            discarded.append(sample)
    # Join tsv read and input list mfr_Collections
    for format_type in mfrs_found_tsv:
        mfrs_found[format_type].mfrs.update(mfrs_found_tsv[format_type].mfrs)

    # Process samples and build
    for id_str in mfr_samples:
        mfr = mfr_samples[id_str]
        print('MultiFileSample: ', id_str, 'format:', mfr.format_name)
        print(mfr.get_files())
    return mfrs_found, discarded


def _determine_mfr_type_from_sample_list(sample_list_file, collections,

                                         format_dict, target_formats,
                                         use_common_fields_as_id):
    mfrs_found = dict()
    sample_list_info = sanitycheck_and_load_sample_tsv(sample_list_file)
    for sample_info in sample_list_info:
        found_format_or_is_leftover = False
        paired = True
        if sample_info['read2']:
            samples = [sample_info['read1'], sample_info['read2']]
        else:
            samples = [sample_info['read1']]
            paired = False
        sample_name = sample_info['sample_name']
        for target_format_type in target_formats:
            if found_format_or_is_leftover:
                break
            mfrs = []
            seq_files = []
            for sample in samples:
                seq_exts = format_dict[target_format_type]['main_exts']
                format_list = format_dict[target_format_type]['formats']
                # Select which mfr type to try
                mfr_type = collections[target_format_type].mfr_type
                try:
                    zip_exts = (
                            format_dict[target_format_type]['secondary_exts'])
                except KeyError:
                    zip_exts = []

                used_ext = test_extension(sample, seq_exts)
                if not used_ext:
                    continue
                    raise UnknownExtensionError(
                            'Extension not recognized\n%s' % sample)
                sample_string, zipped = sample.split(used_ext)
                if zipped and not test_extension(zipped, zip_exts):
                    continue
                path = os.path.dirname(sample_info["read1"])
                sample = os.path.basename(sample_string)
                # Get first section of file.filename for ID
                seq_file = SequenceFile(sample, path=path,
                                        ext=used_ext, ext2=zipped, add_exts=[])
                seq_files.append(seq_file)
                # Check if known format
                mfrs.append(find_format(seq_file, format_list, mfr_type))
            if target_format_type not in mfrs_found:
                mfrs_found[target_format_type] = (
                        collections[target_format_type]())
            if all(issubclass(type(mfr), MultiFileReads) for mfr in mfrs):
                mfr = mfrs[0]
                mfr.set_meta_info(sample_info['meta_info'])
                if sample_name:
                    mfr.set_immutable_id(sample_name)
                elif use_common_fields_as_id:
                    mfr.use_common_fields_as_id()
                if paired:
                    try:
                        mfr.add_files(seq_files[1])
                    except:
                        raise Exception("Cannot yet deal with missmatched"
                                        " formats in samplesheet row")
                mfrs_found[target_format_type].add(mfr, seq_file)
            else:
                sample_split = []
                read_id_1 = os.path.basename(seq_files[0].filename)
                if paired:
                    read_id_2 = os.path.basename(seq_files[1].filename)
                if not sample_name and not paired:
                    sample_name = read_id_1
                elif not sample_name:
                    offset = 0
                    sample_split = ['', '']
                    sample_len = len(read_id_1)
                    if sample_len != len(read_id_2):
                        raise Exception('Error handling fallback determination'
                                        ' of sample_id.\n'
                                        ' Reads have different length\n'
                                        ' %s\n%s' % (read_id_1, read_id_2))

                    for i in range(len(read_id_1)):
                        try:
                            if offset == 2:
                                break
                            if read_id_1[i] == read_id_2[i]:
                                sample_split[offset] += read_id_1[i]
                            else:
                                offset += 1
                        except IndexError as err:
                            print('name:', read_id_1, 'i:', i,
                                  'samp:', read_id_2,
                                  '\n', err, file=sys.stderr)
                    sample_name = ''.join([x.strip('._+-')
                                          for x in sample_split])
                settings_info = list((format_list[0]).values())[0]
                mfr = IlluminaMFR(
                        settings_info['format'].items(),
                        format_name='unknown_fastq')
                mfr.set_immutable_id(sample_name)
                mfr.sample_dict['read1'] = seq_files[0]
                mfr.sample_dict['read2'] = seq_files[1]
                mfr.field_names = ['sample_prefix', 'sample_suffix']
                if sample_split:
                    mfr.non_var_dict = dict(zip(mfr.field_names,
                                                sample_split))
                else:
                    mfr.non_var_dict['sample_name'] = sample_name
                # mfrs_found[target_format_type].leftovers.add(sample, path,
                #                                             zipped, used_ext)
                mfrs_found[target_format_type].add(mfr, seq_file)
            found_format_or_is_leftover = True
        if not found_format_or_is_leftover:
            raise Exception("Error: Filetype in sample list file"
                            " not recognized!")
    return mfrs_found


class MFR_Collection:
    mfrs = {}
    mfr_type = MultiFileReads

    def __init__(self):
        self.mfrs = {}

    def add(self, mfr, seq_file):
        pass

# <<<---- Illumina FastQ Collections ---->>>


class Illumina_MFR_Collection(MFR_Collection):
    mfrs = {}
    leftovers = None
    mfr_type = IlluminaMFR

    def __init__(self):
        super().__init__()
        self.leftovers = Leftovers()
        self.mfrs = {}

    def add(self, mfr, seq_file):
        if mfr.id_string not in self.mfrs:
            self.mfrs[mfr.id_string] = mfr
        else:
            self.mfrs[mfr.id_string].add_files(seq_file)

    def flatten_rename(self, newIDPrefix='S', start_index=1):
        '''
        '''
        result = {}
        index = start_index
        for mfr_id in sorted(self.mfrs.keys()):
            mfr = self.mfrs[mfr_id]
            for samp_dict in mfr.get_files():
                read_num = len([x for x in samp_dict if x[:-1] == 'read'])
                if read_num == 1:
                    # if 'read1' not in samp_dict:
                    #    print(samp_dict, '\n', mfr.format_name)
                    file = samp_dict['read1']
                    sample = SampleInfoSingle(mfr.id_string,
                                              os.path.abspath(str(file)),
                                              [file.filename],
                                              ext=file.ext, zip_ext=file.ext2,
                                              meta_info=mfr.meta_info)
                elif read_num == 2:
                    file1 = samp_dict['read1']
                    file2 = samp_dict['read2']
                    sample = SampleInfoPaired(os.path.abspath(str(file1)),
                                              os.path.abspath(str(file2)),
                                              mfr.id_string,
                                              [file1.filename, file2.filename],
                                              [file1.ext, file2.ext],
                                              [file2.ext2, file2.ext2],
                                              meta_info=mfr.meta_info)
                else:
                    raise AmbigiousPairedReadsError(
                            'To many files map together:\n%s' %
                            repr(mfr).replace('>,', '>,\n'))
                result['%s%i' % (newIDPrefix, index)] = sample
                index += 1
        return result

    def flatten_naive(self):
        result = {}
        for mfr_id in sorted(self.mfrs.keys()):
            mfr = self.mfrs[mfr_id]
            for samp_dict in mfr.get_files():
                id_values = []
                read_num = len([x for x in samp_dict if x[:-1] == 'read'])
                for field in mfr.field_names:
                    if field in mfr.var_fields:
                        id_values.append(samp_dict[field])
                    if field in mfr.common_fields_dict:
                        id_values.append(mfr.common_fields_dict[field])
                    if field == 'read1' and read_num == 1:
                        if mfr.default_field_values['read']:
                            id_values.append(mfr.default_field_values['read'])
                if mfr.immutable_id:
                    id_string = mfr.id_string
                else:
                    id_string = mfr.main_sep.join(id_values)
                add_info = dict((x, y) for x, y in samp_dict.items()
                                if x[:-1] != 'read')
                add_info.update(mfr.common_fields_dict)
                add_info['format'] = mfr.format_name
                if read_num == 1:
                    # if 'read1' not in samp_dict:
                    #    print(samp_dict, '\n', mfr.format_name)
                    file = samp_dict['read1']
                    sample = SampleInfoSingle(id_string,
                                              os.path.abspath(str(file)),
                                              [file.filename],
                                              ext=file.ext, zip_ext=file.ext2,
                                              add_info=add_info,
                                              meta_info=mfr.meta_info)
                elif read_num == 2:
                    file1 = samp_dict['read1']
                    file2 = samp_dict['read2']
                    sample = SampleInfoPaired(os.path.abspath(str(file1)),
                                              os.path.abspath(str(file2)),
                                              id_string,
                                              [file1.filename, file2.filename],
                                              [file1.ext, file2.ext],
                                              [file2.ext2, file2.ext2],
                                              add_info=add_info,
                                              meta_info=mfr.meta_info)
                else:
                    raise AmbigiousPairedReadsError(
                            'To many files map together:\n%s' %
                            repr(mfr).replace('>,', '>,\n'))
                result[id_string] = sample
        return result

    def get_samples(self):
        result = {}
        for mfr_id in sorted(self.mfrs.keys()):
            mfr = self.mfrs[mfr_id]
            container = []
            for samp_dict in mfr.get_files():

                read_num = len([x for x in samp_dict if x[:-1] == 'read'])
                add_info = dict((x, y) for x, y in samp_dict.items()
                                if x[:-1] != 'read')
                add_info.update(mfr.common_fields_dict)
                add_info['format'] = mfr.format_name
                if read_num == 1:
                    # if 'read1' not in samp_dict:
                    #    print(samp_dict, '\n', mfr.format_name)
                    try:
                        file = samp_dict['read1']
                        sample = SampleInfoSingle(mfr.id_string,
                                                  os.path.abspath(str(file)),
                                                  [file.filename],
                                                  ext=file.ext,
                                                  zip_ext=file.ext2,
                                                  add_info=add_info,
                                                  meta_info=mfr.meta_info)
                    except KeyError as err:
                        eprint('get_verbose_samples(): cannot find read1\n '
                               'culprit: %s' % repr(mfr).replace('>,', '>,\n'))

                elif read_num == 2:
                    file1 = samp_dict['read1']
                    file2 = samp_dict['read2']
                    sample = SampleInfoPaired(os.path.abspath(str(file1)),
                                              os.path.abspath(str(file2)),
                                              mfr.id_string,
                                              [file1.filename, file2.filename],
                                              [file1.ext, file2.ext],
                                              [file2.ext2, file2.ext2],
                                              add_info=add_info,
                                              meta_info=mfr.meta_info)
                else:
                    raise AmbigiousPairedReadsError(
                            'To many files map together:\n%s' %
                            repr(mfr).replace('>,', '>,\n'))
                container.append(sample)
            result[mfr_id] = container
        return result


class Leftovers:
    samples = {}
    delims = []
    num_files = 0

    def __init__(self, delims=['_', '.', '+']):
        self.samples = {}
        self.delims = delims
        self.num_files = 0

    def add(self, sample, path, zipped, used_ext):
        for delim in self.delims:
            sample_delim_split = sample.split(delim)
            sample_ID = sample_delim_split[0]
            if len(sample_delim_split) > 1:
                break
        # Sample not seen before
        num = len(self.samples)
        if sample_ID not in self.samples:
            lcp = 0
            self.samples[sample_ID] = ([sample], [path], lcp,
                                       num, [zipped], [used_ext])
            num += 1
        # Sample already seen
        else:
            (prev_sams, prev_paths, prev_lcp, old_num, zippeds, exts) = (
                    self.samples[sample_ID])
            # Use longest common prefix of files to determine read type
            # Note: Not safe if reads are fractured between different flow cell
            # lanes or tiles... this will break
            lcp = len(longest_common_prefix(prev_sams[0], sample))
            self.samples[sample_ID] = (
                    prev_sams+[sample],
                    prev_paths+[path],
                    lcp, old_num, zippeds+[zipped],
                    exts+[used_ext])

    def process_leftovers(self, rename=True, rename_start_index=1):
        results = dict()
        for id_str in self.samples:
            sams, paths, lcp, num, zippeds, exts = self.samples[id_str]
            num += rename_start_index
            identical = all(lcp == len(x) for x in sams)
            if len(sams) == 2 and not identical:
                # print(sams)
                pair = [x[lcp] for x in sams]
                index = pair.index('1')
                ord_ext = [exts[index]]
                ord_zippeds = [zippeds[index]]
                read1 = os.path.join(paths[index], sams[index]
                                     + exts[index]
                                     + zippeds[index] if zippeds else '')
                read1 = os.path.abspath(read1)
                index = pair.index('2')
                ord_ext += [exts[index]]
                ord_zippeds += [zippeds[index]]
                read2 = os.path.join(paths[index],
                                     sams[index]
                                     + exts[index]
                                     + zippeds[index] if zippeds else '')
                read2 = os.path.abspath(read2)
                final_id_string = ('S%i' % num) if rename else id_str
                results[final_id_string] = SampleInfoPaired(
                        read1, read2, id_str,
                        exts=ord_ext,
                        zip_exts=ord_zippeds)
            elif len(sams) == 1 or identical:
                read1 = os.path.join(paths[0],
                                     sams[0]
                                     + exts[0]
                                     + zippeds[0] if zippeds else '')
                final_id_string = ('S%i' % num) if rename else id_str
                results[final_id_string] = SampleInfoSingle(
                        id_str, read1, ext=exts[0],
                        zip_ext=zippeds[0])
            else:
                # Here goes missing logic to deal with flow cell lanes and co
                # print(sams)
                raise AmbigiousPairedReadsError(
                    'Error: Found %i Files for Sample. Expected 1 or 2\n'
                    'Files for id: %s\n%s\n'
                    ' Flow cell Logic is currently missing'
                    '' % (len(sams), id_str, '\n'.join(sams)))
        return results

#  <<<------- Ion-Torrent Bam Collection -------->>>


class IonTorrent_MFR_Collection(MFR_Collection):
    mfrs = {}
    mfr_type = IonTorrentMFR
    leftovers = None

    def __init__(self):
        super().__init__()
        self.leftovers = Leftovers()

#  <<<------- PacBio h5 and meta.xml Collection ------->>>


class PacBio_MFR_Collection(MFR_Collection):
    mfrs = {}
    mfr_type = PacBioMFR
    leftovers = None

    def __init__(self):
        super().__init__()
        self.mfrs = {}

    def add(self, mfr, seq_file):
        if mfr.id_string not in self.mfrs:
            self.mfrs[mfr.id_string] = mfr
        else:
            self.mfrs[mfr.id_string].add_files(seq_file)

    def get_samples(self):
        results = {}
        for mfr_id in sorted(self.mfrs.keys()):
            mfr = self.mfrs[mfr_id]
            # container = []
            for samp_dict in mfr.get_files():
                bax_num = len([x for x in samp_dict if x in '123'])
                file_dict = dict((x if x not in '123'
                                  else 'bax%s' % x, str(y))
                                 for x, y in samp_dict.items()
                                 if x in ['bas', 'metadata']
                                 or y.ext == '.bax' and x in '123')
                add_info = dict((x, str(y)) for x, y in samp_dict.items()
                                if x not in ['1', '2', '3',
                                             'bas', 'metadata'])
                add_info.update(mfr.common_fields_dict)
                add_info['format'] = mfr.format_name
                file_dict['add_info'] = add_info
                if bax_num == 0 and ('bas' in file_dict):
                    # if 'read1' not in samp_dict:
                    #    print(samp_dict, '\n', mfr.format_name)
                    try:
                        file = samp_dict
                        sample = PacBioSampleInfoRS(mfr.id_string,
                                                    **file_dict)
                    except KeyError as err:
                        eprint('get_verbose_samples(): cannot find read1\n '
                               'culprit: %s' % repr(mfr).replace('>,', '>,\n'))
                    results[mfr_id] = sample
                elif bax_num == 3:
                    sample = PacBioSampleInfoRS_II(mfr.id_string,
                                                   **file_dict)
                    results[mfr_id] = sample
                else:
                    # eprint('To many/less files map together:\n%s' %
                    #        repr(mfr).replace('>,', '>,\n'))
                    pass
        return results


def parse_reference_info(reference_list):
    '''Parsing reference file.filenames into reference info objects

    Turns list of input files into reference data containers.
    It determines reference ids, which compression extension is used (if any),
    and which nucleotide sequence file extension is used.

    It only accepts files that end in <seqf_ext>(otional:<comp_ext>)

    Args:
        reference_list(:obj:`list` of :obj:`str`): list of filenames

    Returns:
        :obj:`list` of :obj:`ReferenceInfo`

    Raises:
        UnknownExtensionError: If sequence file extension unknown
    '''
    results = []
    zip_exts = ['.gz', '.xz', '.bz2', '.lzma', '.lzo', '.lz', '.rz']
    seq_exts = ['.fasta', '.fastq', '.fas', '.fna', '.fnq', '.fa']
    for ref, num in zip(reference_list, range(1, len(reference_list)+1)):
        used_ext = test_extension(ref, seq_exts)
        if not used_ext:
            continue
            # raise UnknownExtensionError(
            #     'Extension not recognized\n%s' % ref)
        ref_id, zipped = ref.split(used_ext)
        if zipped and not test_extension(zipped, zip_exts):
            continue
        ref_id = os.path.basename(ref_id)
        results.append(('G%i' % num,
                        ReferenceInfo(ref_id, os.path.abspath(ref),
                                      used_ext, zipped)))
    return results


def find_format(file, formats, mfr_type):
    '''Checks which naming scheme fits
    '''
    out = file
    for format_raw in formats:
        format_name = list(format_raw.keys())[0]
        settings_dict = list(format_raw.values())[0]
        mfr = mfr_type(settings_dict['format'].items(),
                       main_sep=settings_dict['main_sep'],
                       format_name=format_name)
        try:
            mfr.add_files(out)
            id_sep = (settings_dict['id_sep']
                      if 'id_sep' in settings_dict else None)
            id_fields = (settings_dict['id_fields']
                         if 'id_fields' in settings_dict else None)
            if id_sep is not None or id_fields:
                mfr.update_id_field(id_fields=id_fields,
                                    id_sep=id_sep)
        except FormatMismatch:
            # eprint('%s\n is not format: %s' % (file, format_name))
            # eprint(mfr.regex_str)
            continue
        return mfr
    return out


def get_formats_from_file(format_yaml):
    '''Loads a Yaml file containing file naming standards

    Args:
        format_yaml(:obj:`str`): path of format yaml file to be read
    '''
    with open(format_yaml, 'r') as format_fh:
        return ordered_load(format_fh)


def sanitycheck_and_load_sample_tsv(samples_tsv, sep='\t'):
    '''loads sample list from tsv file and sanity checks the tsv

    The file neeeds to adhere to the followoing layout.
    lines begining with a hash symbol should not be included in the tsv
    #<Header all one line broken here for readability>
    SampleName<tab>Read1<tab>Read2<tab>your_meta_info_name1<tab>meta_info_name2
        <tab>and_so_on
    #<Samples>
    CHicG-12<tab>path2read1<tab>path2read2<tab>recovered<tab>group2
    <tab>path2read1<tab>path2read2<tab>not_recovered<tab>group1
    <tab>path2single_end_read<tab><tab>recoverd<tab>group1

    Sample name and READ2 are optional fields
    If read2 is ommited it treats the sample as single end data

    Args:
        samples_tsv(:obj:`str`): path/to/<sample_list>.tsv
    Kwargs:
        sep(:obj:`str` default '\t'[tab]): separator used by file
    Raises:
        SanityCheckFailedError: If sanity check fails
        IOError: If opening of file fails
        csv.Error: If csv library encounters an error

    '''
    output_list = []
    with open(samples_tsv) as sample_tsv_fh:
        error_and_reason = ''
        sample_name_dict = dict()
        has_header = csv.Sniffer().has_header(sample_tsv_fh.read(1024))
        sample_tsv_fh.seek(0)
        tsv_reader = csv.reader(sample_tsv_fh, delimiter=sep)
        if has_header:
            header = next(tsv_reader)

            linenum = 1
        else:
            linenum = 0

        for row in tsv_reader:
            linenum += 1
            line_error = []
            sample_name = row[0]
            # Check uniqueness of sample name if provided
            if sample_name:
                if sample_name in sample_name_dict:
                    line_error.append("sample name not unique: <%s>"
                                      " already used in line %i" % (
                                          sample_name,
                                          sample_name_dict[sample_name]))
                else:
                    sample_name_dict[sample_name] = linenum
            # Check if read files can be found
            read1 = row[1]
            read2 = row[2]
            if not read1:
                line_error.append("read1 field is missing")
            if not os.path.isfile(read1):
                line_error.append("read1: <%s> file can not be found" % read1)
            if read2 and not os.path.isfile(read2):
                line_error.append("read2: <%s>"
                                  " file can not be found" % read2)
            # Construct meta data dictonary
            if has_header:
                meta_dict = dict(zip(header[3:], row[3:]))
            else:
                meta_dict = dict(zip(['metavar_%i' % (i+1)
                                      for i in range(len(row[3:]))],
                                     row[3:]))
            if line_error:
                error_and_reason += '\n    '.join(
                        ['\nError in line %i of parsed file' % linenum]
                        + line_error)
                error_and_reason += '\n'
            elif not error_and_reason:
                info = dict()
                info["sample_name"] = sample_name
                info["read1"] = read1
                info["read2"] = read2
                info["meta_info"] = meta_dict
                output_list.append(info)
    if error_and_reason:
        raise SanityCheckFailedError(error_and_reason)
    return output_list


def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    '''Stolen function to get ordered dictionary from yaml

    https://stackoverflow.com/questions/5121931/
    in-python-how-can-you-load-yaml-mappings-as-ordereddicts/21048064#21048064
    https://stackoverflow.com/users/650222/coldfix
    '''

    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)


def longest_common_prefix(str1, str2):
    '''longest common prefix of two strings
    '''
    return [i[0]
            for i in takewhile(lambda x: (len(set(x)) == 1),
                               zip(str1, str2))]
