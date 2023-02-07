from os import path
import csv

from dataclasses import dataclass
from pybedtools import BedTool

from utils.file_system import write_to_text_file, FolderCreator
from utils.exceptions import OutputError, FolderCreatorError

from primer_designer import PrimerDesigner

@dataclass
class OutputFilesData:
    dir: str
    def fields(self):
        return list(self.__dataclass_fields__.keys())

@dataclass
class SlicerOutputData(OutputFilesData):
    bed: str = ''
    fasta: str = ''

@dataclass
class PrimerOutputData(OutputFilesData):
    bed: str = ''
    csv: str = ''
    
@dataclass
class PrimerDesignerOutputData(OutputFilesData):
    fn: str = ''
    csv: str = ''
    
    
@dataclass
class IpcressOutputData(OutputFilesData):
    stnd: str = ''
    err: str = ''

@dataclass
class DesignOutputData(SlicerOutputData, PrimerOutputData, IpcressOutputData):
    p3_bed: str=''


def timestamped_dir(prefix):
    try:
        FolderCreator.create_timestamped(prefix)
    except FolderCreatorError as err:
        raise OutputError(f'Error creating folder: {err}')
    return FolderCreator.get_dir()

def write_slicer_output(dir_prefix, slices) -> SlicerOutputData:
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


def export_slices_to_csv(slices, dir):
    PRIMER3_OUTPUT_CSV = 'p3_output.csv'

    headers = ['primer', 'sequence', 'chr', 'primer_start', 'primer_end', 'tm', 'gc_percent', 
        'penalty', 'self_any_th', 'self_end_th', 'hairpin_th', 'end_stability']
    rows = construct_csv_format(slices, headers)

    csv_path = path.join(dir, PRIMER3_OUTPUT_CSV)

    with open(csv_path, "w") as p3_fh:
        p3_out = csv.DictWriter(p3_fh, fieldnames=headers)
        p3_out.writeheader()
        p3_out.writerows(rows)

        return csv_path
    
def export_primer_designer_to_csv(primer_designer, output_data: PrimerDesignerOutputData):
    csv_path = path.join(output_data.dir, output_data.fn)
    fields = primer_designer.get_fields()
    with open(csv_path, "w") as f:
        f_out = csv.DictWriter(f,fieldnames=fields)
        f_out.writeheader()
        f_out.writerows(primer_designer.primer_pairs)
        return csv_path

def construct_csv_format(slices, headers):
    rows = []

    for slice_data in slices:
        primers = slice_data['primers']
        for primer in primers:
            primers[primer]['primer'] = primer
            primers[primer]['chr'] = slice_data['chrom']
 
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


def export_to_bed(bed_rows, dir):
    PRIMER_OUTPUT_BED = 'p3_output.bed'

    p3_bed = BedTool(bed_rows)
    bed_path = path.join(dir, PRIMER_OUTPUT_BED)
    p3_bed.saveas(bed_path)

    return bed_path

def write_primer_output(
    prefix = '',
    primers = [],
    existing_dir = '',
    ) -> PrimerOutputData:
    
    if existing_dir:
        dir = existing_dir
    else:
        dir = timestamped_dir(prefix)

    result = PrimerOutputData(dir)

    bed_rows = construct_bed_format(primers)

    result.bed = export_to_bed(bed_rows, dir)
    result.csv = export_slices_to_csv(primers, dir)
    result.dir = dir

    print('Primer files saved:', result.bed, result.csv)

    return result

def write_primer_designer_output(
    prefix = '',
    primer_designer = PrimerDesigner,
    existing_dir = '',
    ) -> PrimerOutputData:
    
    if existing_dir:
        dir = existing_dir
    else:
        dir = timestamped_dir(prefix)

    result = PrimerDesignerOutputData(dir,r'primer_designer.csv')

    result.csv = export_primer_designer_to_csv(primer_designer,result)

    print(f'Primer Designer files saved:{result.csv}')

    return result

def write_ipcress_input(dir, formatted_primers) -> str:
    INPUT_FILE_NAME = 'ipcress_primer_input'

    file_path = write_to_text_file(dir, formatted_primers, INPUT_FILE_NAME)

    return file_path

def write_ipcress_output(stnd = '', err = '', existing_dir = '') -> IpcressOutputData:
    IPCRESS_OUTPUT_TXT = 'ipcress_output'
    
    result = IpcressOutputData(existing_dir)
    
    result.stnd = write_to_text_file(existing_dir, stnd, IPCRESS_OUTPUT_TXT)
    result.err = write_to_text_file(existing_dir, err, IPCRESS_OUTPUT_TXT+"_err")
    
    return result
    
    
 

