from typing import List
import pandas as pd

from ranking.ranking_validator import RankCriteriaValidator, RankingCriteria
from primer.write_primer_output import _get_primers_dataframe


class Ranker:
    def __init__(self, config: dict):
        self._ranking_order: List[RankingCriteria] = RankCriteriaValidator(config).validated_criteria()
    
    def rank(self, primer_type:str, primer_pairs=[]) -> pd.DataFrame:

        if not self._ranking_order:
            return

        primer_pairs_df = _get_primers_dataframe(primer_pairs, primer_type)

        columns_to_sort = [column.column for column in self._ranking_order]

        is_ascending = [column.is_ascending for column in self._ranking_order]

        primer_pairs_df.sort_values(by=columns_to_sort, ascending=is_ascending, inplace=True)

        return primer_pairs_df