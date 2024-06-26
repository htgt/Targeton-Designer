import unittest
import logging
from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase
from textwrap import dedent
from io import StringIO

from utils.exceptions import Primer3Error
from primer.primer_pair import PrimerPair
from primer.primer3_handle_errors import format_no_primer_pairs_message, handle_primer3_errors, _get_primer3_explain


class CapturingStreamHandler(logging.StreamHandler):
  """Custom StreamHandler to capture logs in memory."""
  def __init__(self):
    super().__init__()
    self.buffer = StringIO()
    self.stream = self.buffer


class TestPrimer3ErrorHandling(TestCase):
    
    def setUp(self):
        # Create a custom stream handler to capture logs
        self.handler = CapturingStreamHandler()
        # Get the logger and set the level to capture warnings (adjust if needed)
        logger = logging.getLogger()
        logger.setLevel(logging.WARNING)  # Capture warnings by default
        logger.addHandler(self.handler)

    def tearDown(self):
        # Remove the handler after each test to reset logging
        logger = logging.getLogger()
        logger.removeHandler(self.handler)

    def test_format_no_primer_pairs_message_flag(self):
        # arrange
        stringency = 1
        primer_explain_flag = 1
        designs = {'PRIMER_LEFT_EXPLAIN': 'considered 1469, GC content failed 769, low tm 1, high tm 657, high hairpin stability 2, ok 40', 
                'PRIMER_RIGHT_EXPLAIN': 'considered 1469, GC content failed 235, low tm 1, high tm 1159, ok 74', 
                'PRIMER_PAIR_EXPLAIN': 'considered 2960, unacceptable product size 2960, ok 0', 
                'PRIMER_PAIR_NUM_RETURNED': 0}

        # act
        result = format_no_primer_pairs_message(stringency, primer_explain_flag, designs)

        # assert
        expected_result = 'Stringency level 1 -- PRIMER_LEFT_EXPLAIN: considered 1469, GC content failed 769, low tm 1, high tm 657, high hairpin stability 2, ok 40; PRIMER_RIGHT_EXPLAIN: considered 1469, GC content failed 235, low tm 1, high tm 1159, ok 74; PRIMER_PAIR_EXPLAIN: considered 2960, unacceptable product size 2960, ok 0'

        self.assertEqual(result, expected_result)

    def test_format_no_primer_pairs_message_no_flag(self):
        # arrange
        stringency = 1
        primer_explain_flag = 0
        designs = {'PRIMER_LEFT_EXPLAIN': 'considered 1469, GC content failed 769, low tm 1, high tm 657, high hairpin stability 2, ok 40', 
                'PRIMER_RIGHT_EXPLAIN': 'considered 1469, GC content failed 235, low tm 1, high tm 1159, ok 74', 
                'PRIMER_PAIR_EXPLAIN': 'considered 2960, unacceptable product size 2960, ok 0', 
                'PRIMER_PAIR_NUM_RETURNED': 0}

        # act
        result = format_no_primer_pairs_message(stringency, primer_explain_flag, designs)

        # assert
        expected_result = 'Stringency level 1 -- No primer pairs returned; add PRIMER_EXPLAIN_FLAG == 1 to config file for more details'
        self.assertEqual(result, expected_result)

    def test_handle_primer3_errors_error(self):
        # arrange
        primer_explain = ['Stringency level 1 -- PRIMER_LEFT_EXPLAIN: considered 1469, GC content failed 769, low tm 1, high tm 657, high hairpin stability 2, ok 40; PRIMER_RIGHT_EXPLAIN: considered 1469, GC content failed 235, low tm 1, high tm 1159, ok 74; PRIMER_PAIR_EXPLAIN: considered 2960, unacceptable product size 2960, ok 0']
        primer_pairs_any = 0

        # act
        with self.assertRaises(Primer3Error) as primer_error:
            handle_primer3_errors(primer_explain, primer_pairs_any)

        # assert
        expected_result = dedent("""\
                                 NO PRIMER PAIRS BUILT BY PRIMER3:
                                 Stringency level 1 -- PRIMER_LEFT_EXPLAIN: considered 1469, GC content failed 769, low tm 1, high tm 657, high hairpin stability 2, ok 40; PRIMER_RIGHT_EXPLAIN: considered 1469, GC content failed 235, low tm 1, high tm 1159, ok 74; PRIMER_PAIR_EXPLAIN: considered 2960, unacceptable product size 2960, ok 0"""
        )
        result = str(primer_error.exception)
        self.assertEqual(result, expected_result)

    def test_handle_primer3_errors_warning(self):

        # arrange
        primer_explain = ['Stringency level 1 -- PRIMER_LEFT_EXPLAIN: considered 1469, GC content failed 769, low tm 1, high tm 657, high hairpin stability 2, ok 40; PRIMER_RIGHT_EXPLAIN: considered 1469, GC content failed 235, low tm 1, high tm 1159, ok 74; PRIMER_PAIR_EXPLAIN: considered 2960, unacceptable product size 2960, ok 0']
        primer_pairs_any = 1

        # act
        handle_primer3_errors(primer_explain, primer_pairs_any)

        # assert
        expected_result = dedent("""\
                                 Warning: No primer pairs built by Primer3 with the following stringencies:
                                 Stringency level 1 -- PRIMER_LEFT_EXPLAIN: considered 1469, GC content failed 769, low tm 1, high tm 657, high hairpin stability 2, ok 40; PRIMER_RIGHT_EXPLAIN: considered 1469, GC content failed 235, low tm 1, high tm 1159, ok 74; PRIMER_PAIR_EXPLAIN: considered 2960, unacceptable product size 2960, ok 0"""
        )
        result = self.handler.buffer.getvalue().strip()
        self.assertEqual(result, expected_result)

    if __name__ == '__main__':
        unittest.main()
