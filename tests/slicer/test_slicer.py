import unittest

from unittest.mock import patch

import pybedtools.helpers
from pybedtools import BedTool
from pyfakefs.fake_filesystem_unittest import TestCase

from slicer.slicer import Slicer


class TestSlicer(TestCase):
    def setUp(self):
        self.slicer = Slicer()
        self.setUpPyfakefs()
        self.bed_file_data = 'chr1\t100\t250\texon1\t.\t+'
        self.fasta_file_data = '>region1_1::chr1:5-10(+)\nAGTCT\n>region1_2::chr1:15-20(+)\nATTTT\n'

    def create_test_files(self):
        self.fs.create_file('/test.bed', contents=self.bed_file_data)
        self.fs.create_file('/test.fa', contents=self.fasta_file_data)

    def test_get_slice_data_named_exon_success(self):
        # arrange
        expected = [
            ('chr1', 50, 260, 'exon1_1', '.', '+'),
            ('chr1', 90, 300, 'exon1_2', '.', '+')
        ]
        bed = BedTool(self.bed_file_data, from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 40
        }

        # act
        actual = self.slicer.get_slice_data(bed, params)

        # assert
        self.assertEqual(actual, expected)

    def test_get_slice_data_unnamed_exon_success(self):
        # arrange
        expected = [
            ('chr1', 50, 260, 'region1_1', '.', '-'),
            ('chr1', 90, 300, 'region1_2', '.', '-')
        ]
        bed = BedTool('chr1\t100\t250\t.\t.\t-', from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 40
        }

        # act
        actual = self.slicer.get_slice_data(bed, params)

        # assert
        self.assertEqual(actual, expected)

    def test_decrement_one_based_starts(self):
        # arrange
        input_file = [['1', '200', '300', 'name', '0', '+']]
        expected_row = ['1', '199', '300', 'name', '0', '+']

        # act
        result = self.slicer.decrement_one_based_starts(input_file, [])

        # assert
        self.assertEqual(expected_row, result[0])

    @patch('pybedtools.BedTool.sequence')
    def test_get_seq_throw_bad_error_BEDToolsError(self, sequence_mock):
        # arrange
        sequence_mock.side_effect = pybedtools.helpers.BEDToolsError('throw bad error',
                                                                     'this is an unexpected error')
        input_slice_bed = BedTool([
            ('chr1', 50, 260, 'region1_1', '.', '-'),
            ('chr1', 90, 300, 'region1_2', '.', '-')
        ])
        input_fasta = '/test.fa'
        expected = "\nCommand was:\n\n\t\nCommand was:\n\n\tthrow bad error\n\nError message was:\nthis is an " \
                   "unexpected error\n\nError message was:\nPyBEDTools exited with err type BEDToolsError. " \
                   "Arguments:\n'this is an unexpected error'"

        expected_seq_options = f"fi='{input_fasta}', s=True, name+=True"

        # act
        with self.assertRaises(pybedtools.helpers.BEDToolsError) as exception_context:
            self.slicer.get_seq(input_slice_bed, input_fasta)

        # assert
        self.assertEqual(str(exception_context.exception), expected)
        self.assertTrue(sequence_mock.called)
        self.assertEqual(f"{sequence_mock.call_args}", f"call({expected_seq_options})")

    @patch('pybedtools.BedTool.sequence')
    def test_get_seq_throw_good_error_success(self, sequence_mock):
        # arrange
        sequence_mock.side_effect = pybedtools.helpers.BEDToolsError('throw good error',
                                                                     '*****ERROR: Unrecognized parameter: -name+ *****')
        input_slice_bed = BedTool([
            ('chr1', 50, 260, 'region1_1', '.', '-'),
            ('chr1', 90, 300, 'region1_2', '.', '-')
        ])
        input_fasta = '/test.fa'
        expected = '\nCommand was:\n\n\tthrow good error\n\nError message was:\n*****ERROR: Unrecognized parameter: ' \
                   '-name+ *****'

        expected_seq_options = f"fi='{input_fasta}', s=True, name=True"

        # act
        with self.assertRaises(pybedtools.helpers.BEDToolsError) as exception_context:
            self.slicer.get_seq(input_slice_bed, input_fasta)

        # assert
        self.assertEqual(str(exception_context.exception), expected)
        self.assertEqual(f"{sequence_mock.call_args}", f"call({expected_seq_options})")
        self.assertEqual(sequence_mock.call_count, 2)


if __name__ == '__main__':
    unittest.main()
