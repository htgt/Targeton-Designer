import re
from typing import List
from Bio import SeqIO

from primer.ensembl import get_seq_from_ensembl_by_coords


class SliceData:
    def __init__(self, name: str, start: str, end: str, strand: str, chrom: str, bases: str):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.chrom = chrom
        self.bases = bases
        self.targeton_id = name[0:4]
        self.designs = []

    def __repr__(self):
        return (f'SliceData({self.name}, {self.targeton_id}, {self.start}, {self.end}, {self.strand}, {self.chrom},'
                f' {self.bases})')

    @property
    def p3_input(self):
        return {
            'SEQUENCE_ID': self.name,
            'SEQUENCE_TEMPLATE': self.bases,
        }

    @property
    def surrounding_region(self) -> str:
        surrounding_band = 1000

        return get_seq_from_ensembl_by_coords(
            chromosome=self.chrom,
            start=int(self.start) - surrounding_band,
            end=int(self.end) + surrounding_band
        )


    @staticmethod
    def parse_fasta(fasta: str) -> List['SliceData']:
        with open(fasta) as fasta_data:
            rows = SeqIO.parse(fasta_data, 'fasta')

            slices = []
            for row in rows:
                # Name::Chr:Start-End(Strand)
                # ENSE00000769557_HG8_1::1:42929543-42929753
                match = re.search(
                    r'^(\w+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', row.id)
                if match:
                    parsed_id = _parse_slice(match)

                    slice = SliceData(
                        parsed_id["name"],
                        parsed_id["start"],
                        parsed_id["end"],
                        parsed_id["strand"],
                        parsed_id["chrom"],
                        str(row.seq),
                    )

                    slices.append(slice)

        return slices

    @staticmethod
    def get_first_pre_targeton(fasta: str) -> 'SliceData':
        with open(fasta) as fasta_data:
            rows = SeqIO.parse(fasta_data, 'fasta')

            first_row = next(rows)

            # Name::Chr:Start-End(Strand)
            # ENSE00000769557_HG8_1::1:42929543-42929753
            match = re.search(r'^(\w+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', first_row.id)

            if match:
                return SliceData(
                    name=match.group(1),
                    start=match.group(4),
                    end=match.group(5),
                    strand=match.group(6),
                    chrom=match.group(2),
                    bases=str(first_row.seq),
                )
            else:
                raise ValueError(f"The sequence ID '{first_row.id}' does not match the expected format.")


def _parse_slice(match) -> dict:
    coord_data = {
        'name': match.group(1),
        'start': match.group(4),
        'end': match.group(5),
        'strand': match.group(6),
        'chrom': match.group(2),
    }
    return coord_data
