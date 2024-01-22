import primer3
import re
from collections import defaultdict
from _collections_abc import dict_keys
from typing import Tuple, List

from Bio import SeqIO
from Bio.Seq import Seq

from utils.exceptions import InvalidConfigError
from utils.file_system import parse_json


class Primer3:
    def __init__(self, user_config: str = None) -> None:
        default_config = './src/primer/primer3.config.json'
        self._config = user_config or default_config

    def get_primers(self, fasta: str) -> List[dict]:
        print('Reading Fasta file')
        design_inputs = self.read_input_fasta(fasta)

        print('Designing primers')
        designs = self.primer3_design(design_inputs)

        print('Naming primers')
        slices = self.locate_primers(designs)

        return slices

    def read_input_fasta(self, fasta: str) -> List[dict]:
        with open(fasta) as fasta_data:
            rows = SeqIO.parse(fasta_data, 'fasta')

            slices = []
            for row in rows:
                # Name::Chr:Start-End(Strand)
                # ENSE00000769557_HG8_1::1:42929543-42929753
                match = re.search(r'^(\w+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', row.id)
                if match:
                    slice_data = self.construct_slice_coord_dict(match)
                    p3_input = {
                        'SEQUENCE_ID': slice_data['name'],
                        'SEQUENCE_TEMPLATE': str(row.seq),
                    }
                    slice_data['p3_input'] = p3_input
                    slices.append(slice_data)

        return slices

    def primer3_design(self, primer3_inputs: list) -> List[dict]:
        config_data = self._get_config_data()

        for slice_data in primer3_inputs:
            primer3_input = slice_data['p3_input']

            design = primer3.bindings.designPrimers(primer3_input, config_data)
            slice_data['design'] = design

        return primer3_inputs

    def locate_primers(self, designs: list) -> List[dict]:
        slice_designs = []

        for slice_data in designs:
            design = slice_data['design']
            primer_keys = design.keys()

            primers = self._build_primers_dict(design, primer_keys, slice_data)

            del slice_data['design']
            slice_data['primers'] = primers
            slice_designs.append(slice_data)

        return slice_designs

    def _build_primers_dict(self, design, primer_keys: dict_keys, slice_data: dict) -> defaultdict(dict):
        primers = defaultdict(dict)

        for key in primer_keys:
            primer_details = self.capture_primer_details(key)

            if primer_details:
                libamp_name = self.name_primers(primer_details, slice_data['strand'])
                primer_name = slice_data['name'] + "_" + libamp_name + "_" + primer_details['pair']

                primers[primer_name] = self._build_primer_loci(
                    primers[primer_name], key, design,
                    primer_details, slice_data)

        return primers

    def _build_primer_loci(
            self, primer, key, design, primer_details, slice_data: dict
    ) -> dict:

        primer_field = primer_details['field']

        primer[primer_field] = design[key]
        primer['side'] = primer_details['side']

        if primer_field == 'coords':
            primer_coords = self.calculate_primer_coords(
                primer_details['side'], design[key], slice_data['start'])

            primer['primer_start'] = primer_coords[0]
            primer['primer_end'] = primer_coords[1]
            primer['strand'] = self.determine_primer_strands(
                primer_details['side'], slice_data['strand']
            )

        return primer

    def _get_config_data(self) -> dict:
        try:
            config_data = parse_json(self._config)
        except FileNotFoundError:
            raise
        except Exception:
            raise InvalidConfigError('Primer3 config file is not a correct JSON')

        return config_data

    @staticmethod
    def name_primers(primer_details, strand) -> str:
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

    @staticmethod
    def capture_primer_details(primer_name) -> dict:
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

    @staticmethod
    def construct_slice_coord_dict(match) -> dict:
        coord_data = {
            'name': match.group(1),
            'start': match.group(4),
            'end': match.group(5),
            'strand': match.group(6),
            'chrom': match.group(2),
        }
        return coord_data

    @staticmethod
    def calculate_primer_coords(side, coords, slice_start) -> Tuple[int, int]:
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

    @staticmethod
    def determine_primer_strands(side, slice_strand) -> str:
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

    @staticmethod
    def revcom_reverse_primer(seq, strand) -> Seq:
        seq_obj = Seq(seq)

        if strand == '-':
            seq_obj = seq_obj.reverse_complement()

        return seq_obj
