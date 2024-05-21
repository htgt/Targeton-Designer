import primer3

from typing import List
import os

from primer.slice_data import SliceData
from primer.primer3_prepare_config import prepare_p3_config
from primer.primer_pair import PrimerPair, build_primer_pairs
from primer.primer3_handle_errors import format_no_primer_pairs_message, handle_primer3_errors

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class Primer3:
    def __init__(
            self,
            designer_config: dict,
            p3_config: dict
    ) -> None:

        self._p3_config = p3_config
        self._stringency_vector = designer_config.get('stringency_vector', [""])

    def get_primers(self, fasta: str) -> List[PrimerPair]:
        primer_pairs = []

        logger.info('Reading Fasta file')
        slices = SliceData.parse_fasta(fasta)

        for slice in slices:
            slice_primer_pairs = self._get_primer_pairs(slice)
            primer_pairs.extend(slice_primer_pairs)

        return primer_pairs

    def _get_primer_pairs(self, slice_data: SliceData) -> List[PrimerPair]:
        primer_pairs = []
        primer_explain = []

        for stringency in self._stringency_vector:
            designs = self._get_primer3_designs(slice_data.p3_input, stringency)

            number_pairs = designs['PRIMER_PAIR_NUM_RETURNED']
            primer_explain_flag = self._p3_config['PRIMER_EXPLAIN_FLAG']

            # If Primer3 does not return any primer pairs for this stringency
            if not number_pairs:
                msg = format_no_primer_pairs_message(stringency, primer_explain_flag, designs)
                primer_explain.append(msg)

            else:
                built_primer_pairs = build_primer_pairs(designs, slice_data, stringency)
                primer_pairs.extend(built_primer_pairs)

        # If Primer3 did not return any primer pairs for at least one stringency
        if primer_explain:
            handle_primer3_errors(primer_explain, any(primer_pairs))
        # If Primer3 returns pairs but built_primer_pairs does not
        elif not primer_pairs:
            raise ValueError("No primer pairs returned")
        return primer_pairs

    def _get_primer3_designs(self, slice_info: dict, stringency: int) -> dict:
        self._find_kmer_lists()
        config_data = prepare_p3_config(self._p3_config, stringency)
        return primer3.bindings.design_primers(slice_info, config_data)

    def _find_kmer_lists(self):
        if self._p3_config['PRIMER_MASK_TEMPLATE']:
            kmer_path = self._p3_config['PRIMER_MASK_KMERLIST_PATH']

            if not os.path.isdir(kmer_path):
                msg = f"Missing directory with kmer lists required for masking. Expected path: '{kmer_path}'"
                logger.exception(ValueError(msg))
                raise ValueError(msg)

            else:
                kmer_lists_required = ['homo_sapiens_11.list', 'homo_sapiens_16.list']
                kmer_lists_present = os.listdir(kmer_path)

                kmer_lists_all_expected_present = set(kmer_lists_required).issubset(kmer_lists_present)

                if not kmer_lists_all_expected_present:
                    msg = f"Missing kmer list files required for masking. Expected files: '{kmer_path}homo_sapiens_11.list' and '{kmer_path}homo_sapiens_16.list'"
                    logger.exception(ValueError(msg))
                    raise ValueError(msg)
