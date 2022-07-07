from os import path
import csv

from dataclasses import dataclass
from pybedtools import BedTool

from utils.file_system import FolderCreator
from utils.exceptions import OutputError

@dataclass
class OutputFilesData:
    dir: str

@dataclass
class SlicerOutputData(OutputFilesData):
    bed: str = ''
    fasta: str = ''

def timestamped_dir(prefix):
    try:
        FolderCreator.create_timestamped(prefix)
    except FolderCreatorError as err:
        raise OutputError(f'Error creating folder: {err}')
    return FolderCreator.get_dir()

def write_slicer_output(dir_prefix, slices):
    dir = timestamped_dir(dir_prefix)
    result = SlicerOutputData(dir = dir)

    result.bed = write_slicer_bed_output(dir, slices)
    result.fasta = write_slicer_fasta_output(dir, slices)

    print('Slice files saved: ', result.bed, result.fasta)

    return result

def write_slicer_bed_output(dir, slices):
    BED_OUTPUT = 'slicer_output.bed'

    bed_path = path.join(dir, BED_OUTPUT)
    slices.saveas(bed_path)

    return bed_path


def write_slicer_fasta_output(dir, slices):
    FASTA_OUTPUT = 'slicer_output.fasta'

    fasta_path = path.join(dir, FASTA_OUTPUT)
    slices.save_seqs(fasta_path)

    return fasta_path


def export_to_csv(slices, output_dir):
    headers = ['primer', 'sequence', 'tm', 'gc_percent', 'penalty', 'self_any_th', 'self_end_th', 'hairpin_th',
               'end_stability']
    rows = construct_csv_format(slices, headers)

    path = output_dir + '/p3_output.csv'
    with open(path, "w") as p3_fh:
        p3_out = csv.DictWriter(p3_fh, fieldnames=headers)
        p3_out.writeheader()
        p3_out.writerows(rows)

        return path


def construct_csv_format(slices, headers):
    rows = []

    for slice_data in slices:
        primers = slice_data['primers']
        for primer in primers:
            primers[primer]['primer'] = primer

            del primers[primer]['primer_start']
            del primers[primer]['primer_end']
            del primers[primer]['coords']
            del primers[primer]['side']
            del primers[primer]['strand']

            rows.append(primers[primer])

    return rows


def construct_bed_format(slices):
    rows = []
    for slice_data in slices:
        primers = slice_data['primers']
        for primer in primers:
            primer_data = primers[primer]
            # chr,chrStart,chrEnd,name,score,strand
            # Score unknown until iPCRess
            row = [
                slice_data['chrom'],
                primer_data['primer_start'],
                primer_data['primer_end'],
                primer,
                '0',
                primer_data['strand']
            ]
            rows.append(row)
    return rows


def export_to_bed(bed_rows, output_dir):
    p3_bed = BedTool(bed_rows)
    path = output_dir + '/p3_output.bed'
    p3_bed.saveas(path)

    return path

def write_primer_output(prefix = '', primers = [], existing_dir = ''):
    if existing_dir:
        dir = existing_dir
    else:
        dir = timestamped_dir(prefix)

    bed_rows = construct_bed_format(primers)

    bed_path = export_to_bed(bed_rows, dir)
    csv_path = export_to_csv(primers, dir)

    print('Primer files saved:', bed_path, csv_path)
