from collections import defaultdict
from typing import List

import pandas as pd
from os import path

from designer.output_data_classes import PrimerOutputData
from primer.designed_primer import DesignedPrimer
from primer.primer_pair import PrimerPair
from utils.write_output_files import timestamped_dir, export_to_bed


def write_primer_output(
    prefix='',
    primer_pairs=[],
    existing_dir='',
) -> PrimerOutputData:
    if existing_dir:
        export_dir = existing_dir
    else:
        export_dir = timestamped_dir(prefix)

    result = PrimerOutputData(export_dir)

    primer_rows = construct_primer_rows_bed_format(primer_pairs)
    result.bed = export_to_bed(primer_rows, export_dir)

    result.csv = export_primers_to_csv(primer_pairs, export_dir)
    result.dir = export_dir

    print('Primer files saved:', result.bed, result.csv)

    return result


def export_primers_to_csv(primer_pairs: List[PrimerPair], export_dir: str) -> str:
    PRIMER3_OUTPUT_CSV = 'p3_output.csv'
    primers_csv_output_path = path.join(export_dir, PRIMER3_OUTPUT_CSV)

    primers_dataframe = _get_primers_dataframe(primer_pairs)
    primers_dataframe.to_csv(primers_csv_output_path, index=False)

    return primers_csv_output_path


def _get_primers_dataframe(pairs: List[PrimerPair]) -> pd.DataFrame:
    primers_dict = defaultdict(list)

    for pair in pairs:
        for direction in ['forward', 'reverse']:
            primer = getattr(pair, direction)
            primers_dict['primer'].append(primer.name)
            primers_dict['penalty'].append(primer.penalty)
            primers_dict['stringency'].append(primer.stringency)
            primers_dict['sequence'].append(primer.sequence)
            primers_dict['primer_start'].append(primer.primer_start)
            primers_dict['primer_end'].append(primer.primer_end)
            primers_dict['tm'].append(primer.tm)
            primers_dict['gc_percent'].append(primer.gc_percent)
            primers_dict['self_any_th'].append(primer.self_any_th)
            primers_dict['self_end_th'].append(primer.self_end_th)
            primers_dict['hairpin_th'].append(primer.hairpin_th)
            primers_dict['end_stability'].append(primer.end_stability)

        primers_dict['chromosome'].extend([pair.chromosome] * 2)
        primers_dict['pre_targeton_start'].extend([pair.pre_targeton_start] * 2)
        primers_dict['pre_targeton_end'].extend([pair.pre_targeton_end] * 2)

    return pd.DataFrame(primers_dict)


def construct_primer_rows_bed_format(pairs: List[PrimerPair]) -> list:
    primer_rows = []
    for pair in pairs:
        primer_rows.append(create_bed_row_for_primer(primer=pair.forward, chromosome=pair.chromosome))
        primer_rows.append(create_bed_row_for_primer(primer=pair.reverse, chromosome=pair.chromosome))

    return primer_rows


def create_bed_row_for_primer(primer: DesignedPrimer, chromosome: str) -> list:
    primer_row = [
        chromosome,
        primer.primer_start,
        primer.primer_end,
        primer.name,
        '0',
        primer.strand
    ]

    return primer_row
