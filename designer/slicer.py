#!/usr/bin/env python3
from pybedtools import BedTool
from pybedtools.helpers import BEDToolsError
from os import path
from Bio import SeqIO
import csv
import argparse
import sys
import re


def add_modules_to_sys_path():
    BASE_PATH = path.dirname(path.dirname(path.abspath(__file__)))
    sys.path.append(BASE_PATH)

add_modules_to_sys_path()

from utils.exceptions import FileFormatError
from utils.exceptions import SlicerError
from utils.file_system import FolderCreator, check_file_exists

def validate_files(bed, fasta):
    check_file_exists(bed)
    check_file_exists(fasta)

    validate_bed_format(bed)
    validate_fasta_format(fasta)

    validate_bed_content(bed)

    return


def handle_one_based_input(input_bed):
    adjusted_tsv = []
    with open(input_bed) as file:
        tsv = csv.reader(file, delimiter="\t")
        adjusted_tsv = decrement_one_based_starts(tsv, adjusted_tsv) 
    return adjusted_tsv


def decrement_one_based_starts(tsv, new_tsv):
    #BED is only 0-based on the start thus only need to edit column 1
    for row in tsv:
        row[1] = str(int(row[1]) - 1)
        new_tsv.append(row)
    return new_tsv


def _generate_slice_data(exon, exon_name, params):
    slices = []
    start = exon.start - params['flank_5']
    end = start + params['length']
    count = 1
    while end <= (exon.end + params['flank_3']):
        slice_name = f'{exon_name}_{count}'
        slices.append((exon.chrom, start, end, slice_name,
            exon.score, exon.strand))
        start += params['offset']
        end += params['offset']
        count += 1
    return slices


def get_slice_data(bed, params):
    slices = []
    count = 1
    for exon in bed:
        name = exon.name if exon.name != '.' else f'region{count}'
        slices.extend(_generate_slice_data(exon, name, params))
        count += 1
    return slices


def positive_int(arg):
    if int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0')
    return int(arg)


def len_positive_int(arg):
    if 10000 < int(arg) or int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0 and below 10000')
    return int(arg)

def timestamped_dir(arg):
    try:
        FolderCreator.create_timestamped(arg)
    except FolderCreatorError as err:
        raise SlicerError(f'Error creating folder: {err}')
    return FolderCreator.get_dir()

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Get sequence slices for regions in'
                    ' BED file according to parameters specified')
    parser.add_argument('input_bed',
                        help='BED file containing regions of interest')
    parser.add_argument('input_fasta',
                        help='FASTA file to retrieve sequences from')
    parser.add_argument('-1b',
                        help='Declare if input BED/TSV is one-based',
                        action='store_true')
    parser.add_argument('-f5', '--flank_5',
                        help='how far to extend region at 5\' end (default 50nt)',
                        type=int, default=50)
    parser.add_argument('-f3', '--flank_3',
                        help='how far to extend region at 3\' end (default 50nt)',
                        type=int, default=50)
    parser.add_argument('-l', '--length',
                        help='length of each slice (default 210nt)',
                        type=len_positive_int, default=210)
    parser.add_argument('-o', '--offset',
                        help='offset between each slice (default 5nt)',
                        type=positive_int, default=5)
    parser.add_argument('-d', '--dir',
                        help='output directory name to be timestamped'
                        ' (default \'td_output\')',
                        type=timestamped_dir, default='td_output')
    return parser.parse_args(args)


def get_slices(params):
    validate_files(params['input_bed'], params['input_fasta'])
    
    input_bed = params['input_bed']
    if params['1b']:
        input_bed = handle_one_based_input(params['input_bed'])
    
    bed = BedTool(input_bed)
    slice_bed = BedTool(get_slice_data(bed, params))
    # return named, coords slice sequences on specified strand
    seq_options = { 
        "fi"    : params['input_fasta'],
        "s"     : True,
        "name+" : True
    }
    seq = {}
    try:
        seq = slice_bed.sequence(**seq_options)
    except BEDToolsError as bed_err:
        if not re.search(r'\*{5}ERROR:\ Unrecognized parameter: -name\+\ \*{5}', bed_err.args[1]):
            template = "PyBEDTools exited with err type {0}. Arguments:\n{1!r}"
            message = template.format(type(bed_err).__name__, bed_err.args[1])
            raise BEDToolsError(bed_err, message) 
        del seq_options['name+']
        seq_options['name'] = True
        seq = slice_bed.sequence(**seq_options)

    return seq


def validate_bed_format(bed):
    with open(bed) as file:
        tsv_file = csv.reader(file, delimiter='\t')

        line_num = 1
        for line in tsv_file:
            if len(line) < 6:
                raise FileFormatError(f'Unable to read in BED file correctly. Check file format on line {line_num}.')

            line_num = line_num + 1


def validate_fasta_format(fasta):
    with open(fasta) as handle:
        if not any(SeqIO.parse(handle, "fasta")):
            raise FileFormatError('Unable to read in Fasta file correctly. Check file format.')


def validate_bed_content(bed):
    with open(bed) as file:
        tsv_file = csv.reader(file, delimiter='\t')
        line_num = 1
        for line in tsv_file:

            if not re.search(r'^(?:[Cc][Hh][Rr])?([XxYy]|[1-9]|1\d|2[012])$', line[0]):
                raise ValueError(f'Chromosome format incorrect on line {line_num}: {line[0]}')

            if not re.search(r'^\d+$', line[1]):
                raise ValueError(f'Start coordinate format incorrect on line {line_num}: {line[1]}')

            if not re.search(r'^\d+$', line[2]):
                raise ValueError(f'End coordinate format incorrect on line {line_num}: {line[2]}')

            if int(line[2]) < int(line[1]):
                raise ValueError(f'End coordinate must be greater than start coordinate '
                                 f'on line {line_num}. Start: {line[1]} End: {line[2]}')

            if (int(line[2]) - int(line[1])) > 10000:
                raise ValueError(f'Difference between start coordinate and end coordinate '
                                 f'must be less than 10000. On line {line_num} '
                                 f'Difference: {int(line[2]) - int(line[1])}')

            if not line[3]:
                raise ValueError(f'Error with  name field, if no name is supplied please mark '
                                 f'with a \'.\' on line {line_num}: {line[3]}')

            if not line[4]:
                raise ValueError(f'Error with score field, if no score is supplied please mark '
                                 f'with a \'.\' on line {line_num}: {line[4]}')

            if not re.search(r'^[+-]$', line[5]):
                raise ValueError(f'Strand format incorrect on line {line_num}: {line[5]}')

            line_num = line_num + 1


def main(params):
    try:
        slices = get_slices(params)
        slices.saveas(path.join(params['dir'], 'slices.bed'))
        slices.save_seqs(path.join(params['dir'], 'slices.fa'))
        print('Slice files saved')

    except ValueError as valErr:
        raise SlicerError('Error occurred while checking file content: {0}'.format(valErr))
    except FileFormatError as fileErr:
        raise SlicerError('Error occurred while checking file format: {0}'.format(fileErr))
    except FileNotFoundError as fileErr:
        raise SlicerError('Input file not found: {0}'.format(fileErr))
    except Exception as err:
        raise SlicerError('Unexpected error occurred: {0}'.format(err))


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])
    main(vars(parsed_args))
