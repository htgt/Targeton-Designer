#!/usr/bin/env python3

import sys
import re
from os import path
from warnings import warn
from pybedtools import BedTool
from typing import Tuple

from utils.arguments_parser import ParsedInputArguments
from utils.validate_files import validate_files
from utils.write_output_files import (
    write_slicer_output,
    write_primer_output,
    write_ipcress_input,
    write_ipcress_output,
    write_targeton_csv,
    write_scoring_output,
    DesignOutputData,
    IpcressOutputData,
    PrimerOutputData,
    PrimerDesignerOutputData,
    SlicerOutputData,
    TargetonCSVData,
    ScoringOutputData,
)
from slicer.slicer import Slicer
from primer.primer3 import Primer3
from ipcress.ipcress import Ipcress
from adapters.primer3_to_ipcress import Primer3ToIpcressAdapter
from primer_designer import PrimerDesigner

sys.path.insert(
    0, path.abspath(path.join(path.dirname(__file__), '../sge-primer-scoring/src'))
)
from scoring import Scoring


def version_command():
    python_version = sys.version
    version = '0.0.1'

    print('Targeton Designer version: ', version)
    print('Python version: ', python_version)


def slicer_command(args) -> SlicerOutputData:
    validate_files(bed = args['bed'], fasta = args['fasta'])
    slicer = Slicer()
    slices = slicer.get_slices(args)

    return write_slicer_output(args['dir'], slices)


def primer_command(
        fasta = '',
        prefix = '',
        existing_dir = ''
    ) -> PrimerOutputData:

    validate_files(fasta = fasta)
    p3_class = Primer3()
    primers = p3_class.get_primers(fasta)

    primer_result = write_primer_output(
        primers = primers,
        prefix = prefix,
        existing_dir = existing_dir,
    )
 
    return primer_result

def primer_designer_command(
    design_output_data : DesignOutputData,
    primer_designer = PrimerDesigner(),
    prefix = '',
    existing_dir = ''
    ) -> PrimerDesignerOutputData:

    validate_files(p3_csv = design_output_data.p3_csv, score_tsv = design_output_data.scoring_tsv)
    
    primer_designer.prepare_primer_designer(design_output_data)
    primer_designer_result = primer_designer.write_output(
        prefix = prefix,
        existing_dir = existing_dir,
    )
    return primer_designer_result


def primer_for_ipcress(fasta = '', prefix = '', min = 0, max = 0):
    primer_result = primer_command(fasta = fasta, prefix = prefix)

    adapter = Primer3ToIpcressAdapter()
    adapter.prepare_input(primer_result.csv, min, max, primer_result.dir)

    input_path = write_ipcress_input(primer_result.dir, adapter.formatted_primers)

    result = adapter.formatted_primers

    return result


def ipcress_command(params, csv='', existing_dir='') -> IpcressOutputData:
    ipcress_params = params.copy()

    if csv:
        ipcress_params['p3_csv'] = csv
    if existing_dir:
        ipcress_params['dir'] = existing_dir

    validate_files(
        fasta = ipcress_params['fasta'],
        txt = ipcress_params['primers'],
        p3_csv = ipcress_params['p3_csv']
    )

    ipcress = Ipcress(ipcress_params)
    ipcress_result = ipcress.run()

    result = write_ipcress_output(
        stnd = ipcress_result.stnd,
        err = ipcress_result.err,
        existing_dir = ipcress_params['dir']
    )
    result.input_file = ipcress_result.input_file

    return result


def generate_targeton_csv(ipcress_input, bed, dirname, dir_timestamped=False) -> TargetonCSVData:
    bed = BedTool(bed)
    csv_rows = []
    with open(ipcress_input) as fh:
        ipcress_input_data = fh.read()
    for region in bed:
        # corresponding primer pair names will be prefixed by region name
        for primer_pair in re.finditer(rf'^{region.name}\S*', ipcress_input_data, re.MULTILINE):
            csv_rows.append([primer_pair.group(), region.name])
    return write_targeton_csv(csv_rows, dirname, dir_timestamped)


def scoring_command(ipcress_output, mismatch, output_tsv, targeton_csv=None) -> ScoringOutputData:
    scoring = Scoring(ipcress_output, mismatch, targeton_csv)
    scoring.add_scores_to_df()

    result = write_scoring_output(scoring, output_tsv)

    return result


def design_command(args) -> DesignOutputData:
    slicer_result = slicer_command(args)
    primer_result = primer_command(fasta = slicer_result.fasta, existing_dir = slicer_result.dir)
    ipcress_result = ipcress_command(args, csv = primer_result.csv, existing_dir = slicer_result.dir)
    targeton_result = generate_targeton_csv(
        ipcress_result.input_file, args['bed'], slicer_result.dir, dir_timestamped=True
    )
    scoring_output_path = path.join(slicer_result.dir, 'scoring_output.tsv')
    scoring_result = scoring_command(
        ipcress_result.stnd, args['mismatch'], scoring_output_path, targeton_result.csv
    )
    design_result = DesignOutputData(slicer_result.dir)
    # Slicer
    design_result.slice_bed = slicer_result.bed
    design_result.slice_fasta = slicer_result.fasta
    # Primer
    design_result.p3_bed = primer_result.bed
    design_result.p3_csv = primer_result.csv
    # iPCRess
    design_result.ipcress_input = ipcress_result.input_file
    design_result.ipcress_output = ipcress_result.stnd
    design_result.ipcress_err = ipcress_result.err
    # Targeton CSV
    design_result.targeton_csv = targeton_result.csv
    # Scoring
    design_result.scoring_tsv = scoring_result.tsv
    # Primer Designer
    primer_designer_result = primer_designer_command(design_result, existing_dir = slicer_result.dir)   
    design_result.pd_json = primer_designer_result.json
    design_result.pd_csv = primer_designer_result.csv
    
    return design_result


def resolve_command(args):
    command = args['command']

    if command == 'version':
        version_command()
    else:
        if command == 'slicer':
            slicer_command(args)

        if command == 'primer':
            primer_command(fasta = args['fasta'], prefix = args['dir'], args = args)

        if command == 'primer_for_ipcress':
            primer_for_ipcress(fasta = args['fasta'], prefix = args['dir'], min = args['min'], max = args['max'])
            
        if command == 'primer_designer':
            design_output_data = DesignOutputData()
            design_output_data.p3_csv =  args['p3_csv'] 
            design_output_data.scoring_tsv = args['score_tsv']
            primer_designer_command(design_output_data, prefix = args['dir'])

        if command == 'ipcress':
            ipcress_command(args)

        if command == 'generate_targeton_csv':
            generate_targeton_csv(args['primers'], args['bed'], args['dir'])

        if command == 'scoring':
            scoring_command(
                args['ipcress_file'],
                args['scoring_mismatch'],
                args['output_tsv'],
                args['targeton_csv'],
            )

        if command == 'design':
            design_command(args)


def main():
    parsed_input = ParsedInputArguments()
    args = parsed_input.get_args()

    resolve_command(args)


if __name__ == '__main__':
    main()
