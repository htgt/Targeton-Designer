from collections import defaultdict
from typing import List, Tuple

from primer.filter.filter import Filter
from primer.filter.filter_response import FilterResponse, PrimerPairDiscarded
from primer.primer_pair import PrimerPair


class DuplicatesFilter(Filter):
    key: str = 'duplicates'
    value_type: type = bool
    reason_discarded: str = "has duplicated with a higher stringency"

    def apply(self, pairs: List[PrimerPair]) -> FilterResponse:
        pairs_to_keep = []
        pairs_to_discard = []

        pairs_duplicates_grouped = _group_duplicates_pairs(pairs)

        for duplicate_group in pairs_duplicates_grouped:
            pair_max_stringency, others = _take_pair_with_max_stringency_and_others(duplicate_group)
            primer_pairs_to_discard = [PrimerPairDiscarded(pair, DuplicatesFilter.reason_discarded) for pair in others]

            pairs_to_keep.append(pair_max_stringency)
            pairs_to_discard.extend(primer_pairs_to_discard)

        return FilterResponse(pairs_to_keep, pairs_to_discard)


def _take_pair_with_max_stringency_and_others(pairs: List[PrimerPair]) -> Tuple[PrimerPair, List[PrimerPair]]:
    pair_max_stringency = _get_max_stringency_pair(pairs)
    other_pairs = [pair for pair in pairs if pair != pair_max_stringency]

    return pair_max_stringency, other_pairs


def _get_max_stringency_pair(duplicate_group: List[PrimerPair]) -> PrimerPair:
    return max(duplicate_group, key=lambda pair: pair.forward.stringency)


def _group_duplicates_pairs(pairs: List[PrimerPair]) -> List[List[PrimerPair]]:
    groups_duplicates = defaultdict(list)
    for pair in pairs:
        groups_duplicates[pair].append(pair)
    return list(groups_duplicates.values())
