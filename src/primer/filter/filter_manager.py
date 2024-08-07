import sys
from typing import List

from primer.filter.filter import Filter
from primer.filter.hap1_variant_filter import HAP1VariantFilter
from primer.filter.duplicates_filter import DuplicatesFilter
from primer.filter.filter_response import FilterResponse
from primer.primer_pair import PrimerPair

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class FilterManager:

    def __init__(self, apply_filters: dict):

        self.filters: List[Filter] = [DuplicatesFilter(), HAP1VariantFilter()]

        self._filters_to_apply: List[Filter] = []

        if len(apply_filters) == 1 and apply_filters.get("duplicates"):
            logger.info("Using duplicates filter as a default filter.")

        correct_filter_names = [filter_object.key for filter_object in self.filters]
        incorrect_filter_names = set(apply_filters.keys()) - set(correct_filter_names)

        if incorrect_filter_names:
            [logger.info(f"Incorrect filter name {incorrect_name}") for incorrect_name in incorrect_filter_names]
            logger.info("Please check and re-run the command again with the correct filter name in the configuration file.")

            # Terminate the running script until incorrect filter names gets corrected.
            sys.exit(1)

        incorrect_values = [key for key, value in apply_filters.items() if not isinstance(value, bool)]
        if incorrect_values:
            msg: str = f"Wrong value(s) provided for '{', '.join(incorrect_values)}' in config " \
                       "file (only takes true or false). Unable to apply filtering - Exiting programme"
            logger.error(msg)
            sys.exit(1)

        # Filter out filters which are not in apply_filter list.
        self.filters = filter(lambda f: f.key in list(apply_filters.keys()), self.filters)

        for filter_object in self.filters:
            if apply_filters[filter_object.key]:
                self._filters_to_apply.append(filter_object)


    def apply_filters(self, primer_pairs_data: List[PrimerPair]) -> FilterResponse:

        pairs_to_keep = primer_pairs_data
        pairs_to_discard = []

        for _filter in self._filters_to_apply:

            logger.info(f"Filter {_filter.key} is applied.")

            if _filter.key == "HAP1_variant":
                logger.info("Requesting HAP1 variant Web service...")

            filter_response = _filter.apply(pairs_to_keep)

            pairs_to_keep = filter_response.primer_pairs_to_keep
            pairs_to_discard.extend(filter_response.primer_pairs_to_discard)

        if not pairs_to_keep:
            logger.warning("All primer pairs discarded during filtering.")

        return FilterResponse(primer_pairs_to_keep=pairs_to_keep,
                              primer_pairs_to_discard=pairs_to_discard)
