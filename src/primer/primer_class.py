from dataclasses import dataclass
from typing import Tuple, List
from collections import defaultdict
from _collections_abc import dict_keys
from Bio.Seq import Seq
import re

from primer.slice_data import SliceData


@dataclass
class Primer:
    name: str
    bases: str
    start: str
    end: str
    chromosome: str
    penalty: str
    stringency: str
    tm: str
    gc_percent: str
    self_any_th: str
    self_end_th: str
    hairpin_th: str
    end_stability: str


def parse_designs_to_primers(slices: List[SliceData]) -> List[Primer]:
    slice_designs = []

    for slice_data in slices:
        slice_data.primers = {}
        for design in slice_data.designs:

            primer_keys = design.keys()

            primers = build_primers_dict(
                design,
                primer_keys,
                slice_data,
                design['stringency']
            )

            for primer in primers:
                slice_data.primers[primer] = primers[primer]
                slice_designs.append(slice_data)

        del slice_data.designs

    return slice_designs


def build_primer_loci(
        primer,
        key,
        design,
        primer_details,
        slice_data: SliceData,
        primer_name: str,
        primer_pair_id: str,
        stringency: str = "",
) -> dict:
    
    primer_field = primer_details['field']

    primer['primer'] = primer_name
    primer[primer_field] = design[key]

    primer['side'] = primer_details['side']

    if stringency != "":
        primer['stringency'] = stringency
    primer['pair_id'] = primer_pair_id

    if primer_field == 'coords':
        primer_coords = calculate_primer_coords(primer_details['side'],
            design[key], slice_data.start)

        primer['primer_start'] = primer_coords[0]
        primer['primer_end'] = primer_coords[1]
        primer['strand'] = determine_primer_strands(
            primer_details['side'], slice_data.strand)

    return primer

def name_primers(primer_details: dict, strand: str) -> str:
        fwd_primers = {
            'left': 'LibAmpF',
            'right': 'LibAmpR',
        }
        rev_primers = {
            'left': 'LibAmpR',
            'right': 'LibAmpF',
        }
        names = {
            '+': fwd_primers,
            '-': rev_primers,
        }

        primer_name = names[strand][primer_details['side']]

        return primer_name

def capture_primer_details(primer_name: str) -> dict:
        match = re.search(r'^(primer_(left|right)_(\d+))(\_(\S+))?$', primer_name.lower())
        result = {}
        if match:
            primer_id = match.group(1)
            primer_side = match.group(2)
            pair_number = match.group(3)
            primer_field = match.group(5)
            if primer_field is None:
                primer_field = 'coords'
            result = {
                'id': primer_id,
                'side': primer_side,
                'field': primer_field,
                'pair': pair_number
            }

        return result

def calculate_primer_coords(side: str, coords: list, slice_start: str) -> Tuple[int, int]:
        slice_start = int(slice_start)
        left_flank = {
            'start': slice_start,
            'end': slice_start + int(coords[1])
        }

        slice_end = slice_start + int(coords[0])
        right_flank = {
            'start': 1 + slice_end - int(coords[1]),
            'end': 1 + slice_end,
        }

        slice_coords = {
            'left': left_flank,
            'right': right_flank
        }

        start = slice_coords[side]['start']
        end = slice_coords[side]['end']

        return start, end

def determine_primer_strands(side: str, slice_strand: str) -> str:
        positive = {
            'left': '+',
            'right': '-',
        }

        negative = {
            'left': '-',
            'right': '+',
        }

        strands = {
            '+': positive,
            '-': negative,
        }

        return strands[slice_strand][side]

def revcom_reverse_primer(seq: str, strand: str) -> Seq:
        seq_obj = Seq(seq)

        if strand == '-':
            seq_obj = seq_obj.reverse_complement()

        return seq_obj

def build_primers_dict(
        design,
        primer_keys: dict_keys,
        slice_data: dict,
        stringency: str = "",
) -> defaultdict(dict):

        primers = defaultdict(dict)

        for key in primer_keys:
            primer_details = capture_primer_details(key)

            if primer_details:
                libamp_name = name_primers(primer_details, slice_data.strand)
                primer_name = slice_data.name + "_" + libamp_name + "_" + \
                              primer_details['pair']

                primer_name_with_stringency = primer_name + "_str" + stringency.replace(
                    ".", "_")
                primer_pair_id = slice_data.name + "_" + primer_details[
                    'pair'] + "_str" + stringency.replace(".", "_")

                primers[primer_name_with_stringency] = \
                    build_primer_loci(
                        primers[primer_name_with_stringency],
                        key,
                        design,
                        primer_details,
                        slice_data,
                        primer_name,
                        primer_pair_id,
                        stringency,
                    )

        return primers





    