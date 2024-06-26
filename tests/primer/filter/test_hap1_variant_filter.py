import unittest

from primer.designed_primer import DesignedPrimer, Interval
from primer.primer_pair import PrimerPair
from primer.filter.hap1_variant_filter import HAP1VariantFilter
from primer.primer_pair_discarded import PrimerPairDiscarded


class TestHAP1VariantFilter(unittest.TestCase):

    def setUp(self) -> None:
        # There is a HAP1 Variant at 11542 in chromosome 1
        self.primer_with_variant = DesignedPrimer(
            name="primer_with_variant",
            penalty=0.5,
            pair_id="pair_id",
            sequence="ATCGATCG",
            coords=Interval(start=199, end=18),
            primer_start=11540,
            primer_end=11545,
            strand="+",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        self.primer_with_no_variant = DesignedPrimer(
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
        pair_with_variant = PrimerPair(
            pair_id="pair_with_hap1_variant",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=200,
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid")
        pair_with_variant.forward = self.primer_with_variant
        pair_with_variant.reverse = self.primer_with_no_variant

        pair_with_no_variant = PrimerPair(
            pair_id="pair_with_no_variant",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=200,
            stringency=1,
            targeton_id="targeton_id",
            uid="uid")
        pair_with_no_variant.forward = self.primer_with_no_variant
        pair_with_no_variant.reverse = self.primer_with_no_variant

        # Act
        pairs_to_filter = [pair_with_variant, pair_with_no_variant]
        filter_response = HAP1VariantFilter().apply(pairs_to_filter)

        # Assertion
        self.assertEqual(len(filter_response.primer_pairs_to_keep), 1)
        self.assertIn(pair_with_no_variant, filter_response.primer_pairs_to_keep)

        self.assertEqual(len(filter_response.primer_pairs_to_discard), 1)
        self.assertIn(PrimerPairDiscarded(pair_with_variant, reason_discarded=HAP1VariantFilter.reason_discarded),
                      filter_response.primer_pairs_to_discard)
