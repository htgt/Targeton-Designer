import base64
from typing import Tuple, List, Optional
from collections import defaultdict
import re

from primer.filter.hap1 import contain_variant
from primer.slice_data import SliceData


class PrimerPair:
    def __init__(self, pair_id: str, chromosome: str,
                       pre_targeton_start: str,
                       pre_targeton_end: str):
        self.id = pair_id
        self.chromosome = chromosome
        self.pre_targeton_start = pre_targeton_start
        self.pre_targeton_end = pre_targeton_end
        self.forward = {}
        self.reverse = {}

    def __repr__(self):
        return (f"PrimerPair(pair_id='{self.id}', chromosome='{self.chromosome}', "
                f"pre_targeton_start='{self.pre_targeton_start}', "
                f"pre_targeton_end='{self.pre_targeton_end}', "
                f"forward={self.forward}, reverse={self.reverse})")

    def __eq__(self, other):
        if isinstance(other, PrimerPair):
            return (
                    self.chromosome == other.chromosome and
                    self.forward == other.forward and
                    self.reverse == other.reverse
            )
        return False

    def __hash__(self):
        return hash((self.chromosome, self.forward, self.reverse))

    @property
    def contain_hap_one_variant(self) -> bool:
        forward_start, forward_end = self.forward['primer_start'], self.forward['primer_end']
        reverse_start, reverse_end = self.reverse['primer_start'], self.reverse['primer_end']

        return (contain_variant(self.chromosome, forward_start, forward_end) or
                contain_variant(self.chromosome, reverse_start, reverse_end))


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


def name_primers(side: str, strand: str) -> str:
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

    primer_name = names[strand][side]

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


def build_primer_pairs(
        design,
        slice_data: SliceData,
        stringency: str = "",
) -> List[PrimerPair]:
    primer_pairs = []
    primers = defaultdict(dict)

    for key in design.keys():
        primer_details = capture_primer_details(key)

        if primer_details:
            libamp_name = name_primers(primer_details['side'], slice_data.strand)
            primer_name = slice_data.name + "_" + libamp_name + "_" + \
                          primer_details['pair']

            primer_name_with_stringency = primer_name + "_str" + stringency.replace(
                ".", "_")

            make_id = str(slice_data.name + "_" + primer_details['pair'] + "_str" + stringency.replace(".", "_"))

            # Following line will generate unique id for each primer pair based on slice_data, primer_details, and stringency.
            
            primer_pair_id = base64.b64encode(make_id.encode('utf-8')).decode()

            primer = build_primer_loci(
                primers[primer_name_with_stringency],
                key,
                design,
                primer_details,
                slice_data,
                primer_name,
                primer_pair_id,
                stringency,
            )

            pair = _find_pair_by_id(primer_pairs, primer_pair_id)

            if pair is None:

                pair = PrimerPair(primer_pair_id, slice_data.chrom,
                                  slice_data.start, slice_data.end)
                
                primer_pairs.append(pair)

            if libamp_name == "LibAmpF":
                pair.forward = primer
            if libamp_name == "LibAmpR":
                pair.reverse = primer

    return primer_pairs


def _find_pair_by_id(pairs: List[PrimerPair], pair_id: str) -> Optional[PrimerPair]:
    for pair in pairs:
        if pair.id == pair_id:
            return pair
    return None
