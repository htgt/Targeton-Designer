from unittest import TestCase

from primer.designed_primer import DesignedPrimer, Interval
from primer.filter.duplicates_filter import DuplicatesFilter
from primer.primer_pair import PrimerPair
from primer.primer_pair_discarded import PrimerPairDiscarded


class TestDuplicatesFilter(TestCase):

    def setUp(self) -> None:
        self.designed_primer = DesignedPrimer(
            name="primer_with_no_variant",
            penalty=0.5,
            pair_id="pair_id",
            sequence="ATCGATCG",
            coords=Interval(start=199, end=18),
            primer_start=10,
            primer_end=20,
            strand="+",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

    def test_apply_filters(self):
        # Arrange
        pair_max_stringency = PrimerPair(
            pair_id="pair_max_stringency",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=200,
            stringency=1,
            targeton_id="targeton_id",
            uid="uid")
        pair_max_stringency.forward = self.designed_primer
        pair_max_stringency.reverse = self.designed_primer

        pair_min_stringency = PrimerPair(
            pair_id="pair_min_stringency",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=200,
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid")
        pair_min_stringency.forward = self.designed_primer
        pair_min_stringency.reverse = self.designed_primer

        # Act
        pairs_to_filter = [pair_max_stringency, pair_min_stringency]
        filter_response = DuplicatesFilter().apply(pairs_to_filter)

        # Assertion
        self.assertEqual(len(filter_response.primer_pairs_to_keep), 1)
        self.assertTrue(pair_max_stringency, filter_response.primer_pairs_to_keep)

        self.assertEqual(len(filter_response.primer_pairs_to_discard), 1)
        self.assertIn(PrimerPairDiscarded(pair_min_stringency, DuplicatesFilter.reason_discarded),
                      filter_response.primer_pairs_to_discard)
