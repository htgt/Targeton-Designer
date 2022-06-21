from __future__ import with_statement
import unittest
from wsgiref import validate
import pybedtools
import argparse

from pyfakefs.fake_filesystem_unittest import TestCase
from io import StringIO

from designer.slicer import (FileFormatError, len_positive_int, check_file_exists, get_slice_data,
                             positive_int, parse_args, get_slices, validate_bed_content, validate_bed_format,
                             validate_fasta_format)


class TestSlicer(TestCase):

    def setUp(self):
        self.setUpPyfakefs()
        self.bed_file_data = 'chr1\t100\t250\texon1\t.\t+'
        self.fasta_file_data = '>region1_1::chr1:5-10(+)\nAGTCT\n>region1_2::chr1:15-20(+)\nATTTT\n'

    def create_test_files(self):
        self.fs.create_file('/test.bed', contents=self.bed_file_data)
        self.fs.create_file('/test.fa', contents=self.fasta_file_data)

    # def test_validate_files(self):
    #   bed = pybedtools.BedTool('chr1\t100\t250\t.\t.\t+',
    #      from_string=True)
    # fasta = pybedtools.example_filename('test.fa')
    # validate_files(bed, fasta)

    def test_check_file_exists_valid_bed_arg_success(self):
        # arrange
        test_arg = '/test.bed'
        self.create_test_files()

        # act
        check_file_exists(test_arg)

    def test_check_file_exists_valid_fasta_arg_success(self):
        # arrange
        test_arg = '/test.fa'
        self.create_test_files()

        # act
        check_file_exists(test_arg)

    def test_check_file_exists_invalid_args_fail(self):
        # arrange
        test_arg = '/test.bed'
        expected = "Unable to find file: /test.bed"

        # act
        with self.assertRaises(FileNotFoundError) as exception_context:
            check_file_exists(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_format_valid_bed_success(self):
        # arrange
        test_arg = '/test.bed'
        self.create_test_files()

        # act
        validate_bed_format(test_arg)

    def test_validate_bed_format_valid_long_bed_success(self):
        # arrange
        test_arg = '/longtest.bed'
        self.fs.create_file('/longtest.bed', contents='chr1\t100\t250\texon1\t.\t+\tfoo\tbar')

        # act
        validate_bed_format(test_arg)

    def test_validate_bed_format_invalid_csv_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1,100,250,exon1,.,+')
        expected = 'Unable to read in BED file correctly. Check file format on line 1.'

        # act
        with self.assertRaises(FileFormatError) as exception_context:
            validate_bed_format(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_format_invalid_len_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\texon1\t.')
        expected = 'Unable to read in BED file correctly. Check file format on line 1.'

        # act
        with self.assertRaises(FileFormatError) as exception_context:
            validate_bed_format(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_fasta_format_valid_fasta_success(self):
        # arrange
        test_arg = '/test.fa'
        self.create_test_files()

        # act
        validate_fasta_format(test_arg)

    def test_validate_fasta_format_invalid_bed_fail(self):
        # arrange
        test_arg = '/test.bed'
        self.create_test_files()
        expected = 'Unable to read in FastA file correctly. Check file format.'

        # act
        with self.assertRaises(FileFormatError) as exception_context:
            validate_fasta_format(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_valid_bed_success(self):
        # arrange
        test_arg = '/test.bed'
        self.create_test_files()

        # act
        validate_bed_content(test_arg)

    def test_validate_bed_content_valid_bed_alt_success(self):
        # arrange
        test_arg = '/goodtest.bed'
        self.fs.create_file('/goodtest.bed', contents='1\t100\t250\texon1\t.\t+')

        # act
        validate_bed_content(test_arg)

    def test_validate_bed_content_invalid_chr_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='NOCHR\t100\t250\texon1\t.\t+')
        expected = 'Chromosome format incorrect on line 1: NOCHR'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_chr_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='\t100\t250\texon1\t.\t+')
        expected = 'Chromosome format incorrect on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_alpha_start_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\ta\t250\texon1\t.\t+')
        expected = 'Start coordinate format incorrect on line 1: a'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_start_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t\t250\texon1\t.\t+')
        expected = 'Start coordinate format incorrect on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_negative_start_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t-100\t250\texon1\t.\t+')
        expected = 'Start coordinate format incorrect on line 1: -100'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_alpha_end_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\ta\texon1\t.\t+')
        expected = 'End coordinate format incorrect on line 1: a'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_end_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t\texon1\t.\t+')
        expected = 'End coordinate format incorrect on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_negative_end_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t-250\texon1\t.\t+')
        expected = 'End coordinate format incorrect on line 1: -250'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_larger_start_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t1000\t250\texon1\t.\t+')
        expected = 'End coordinate must be greater than start coordinate on line 1. Start: 1000 End: 250'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_big_diff_cood_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t20100\texon1\t.\t+')
        expected = 'Difference between start coordinate and end coordinate must be less than 10000. On line 1 Difference: 20000'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_name_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\t\t.\t+')
        expected = 'Error with name field, if no name is supplied please mark with a \'.\' on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_score_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\texon1\t\t+')
        expected = 'Error with score field, if no score is supplied please mark with a \'.\' on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_invalid_strand_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\texon1\t.\ta')
        expected = 'Strand format incorrect on line 1: a'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_both_strand_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\texon1\t.\t-+')
        expected = 'Strand format incorrect on line 1: -+'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_get_slice_data_named_exon_success(self):
        # arrange
        expected = [
            ('chr1', 50, 260, 'exon1_1', '.', '+'),
            ('chr1', 90, 300, 'exon1_2', '.', '+')
        ]
        bed = pybedtools.BedTool(self.bed_file_data, from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 40
        }

        # act
        actual = get_slice_data(bed, params)

        # assert
        self.assertEqual(actual, expected)

    def test_get_slice_data_unnamed_exon_success(self):
        # arrange
        expected = [
            ('chr1', 50, 260, 'region1_1', '.', '-'),
            ('chr1', 90, 300, 'region1_2', '.', '-')
        ]
        bed = pybedtools.BedTool('chr1\t100\t250\t.\t.\t-', from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 40
        }

        # act
        actual = get_slice_data(bed, params)

        # assert
        self.assertEqual(actual, expected)

    def test_positive_int_positive_arg_success(self):
        # arrange
        test_arg = 1
        expected = 1

        # act
        actual = positive_int(test_arg)

        # assert
        self.assertEqual(actual, expected)

    def test_positive_int_negative_arg_fail(self):
        # arrange
        test_arg = -1
        expected = "Parameter must be above 0"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_positive_int_zero_fail(self):
        # arrange
        test_arg = 0
        expected = 'Parameter must be above 0'

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_len_positive_int_positive_arg_success(self):
        # arrange
        test_arg = 1
        expected = 1

        # act
        actual = len_positive_int(test_arg)

        # assert
        self.assertEqual(actual, expected)

    def test_len_positive_int_zero_fail(self):
        # arrange
        test_arg = 0
        expected = 'Parameter must be above 0 and below 10000'

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            len_positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_len_positive_int_negative_arg_fail(self):
        # arrange
        test_arg = -1
        expected = "Parameter must be above 0 and below 10000"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            len_positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_len_positive_int_large_arg_fail(self):
        # arrange
        test_arg = 10001
        expected = "Parameter must be above 0 and below 10000"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            len_positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_parse_args_valid_args_success(self):
        # arrange
        test_args = ['bed', 'fasta', '-f5', '50',
                     '--length', '200', '--output_bed', 'slices.bed']
        expected = {
            'bed': 'bed',
            'fasta': 'fasta',
            'flank5': 50,
            'flank3': 50,
            'length': 200,
            'offset': 5,
            'output_bed': 'slices.bed'
            }

        # act
        actual = parse_args(test_args)

        # assert
        self.assertEqual(actual.input_bed, expected['bed'])
        self.assertEqual(actual.input_fasta, expected['fasta'])
        self.assertEqual(actual.flank_5, expected['flank5'])
        self.assertEqual(actual.flank_3, expected['flank3'])
        self.assertEqual(actual.length, expected['length'])
        self.assertEqual(actual.offset, expected['offset'])
        self.assertEqual(actual.output_bed, expected['output_bed'])
        self.assertIsNone(actual.output_fasta)

    def test_get_slices(self):
        expected_bed = (
            'chr1\t5\t10\tregion1_1\t.\t+\n'
            'chr1\t15\t20\tregion1_2\t.\t+\n'
        )
        expected_fasta = (
            '>region1_1::chr1:5-10(+)\n'
            'AGTCT\n'
            '>region1_2::chr1:15-20(+)\n'
            'ATTTT\n'
        )
        in_bed = StringIO('chr1\t5\t20\t.\t.\t+')
        in_fasta = pybedtools.example_filename('test.fa')
        params = {
            'input_bed': in_bed,
            'input_fasta': in_fasta,
            'flank_5': 0,
            'flank_3': 0,
            'length': 5,
            'offset': 10
        }
        slices = get_slices(params)
        self.assertEqual(expected_bed, slices.head(as_string=True))
        self.assertEqual(expected_fasta, slices.print_sequence())
        expected_bed = (
            'chr1\t5\t10\texon1_1\t.\t-\n'
            'chr1\t15\t20\texon1_2\t.\t-\n'
        )
        expected_fasta = (
            '>exon1_1::chr1:5-10(-)\n'
            'AGACT\n'
            '>exon1_2::chr1:15-20(-)\n'
            'AAAAT\n'
        )
        params['input_bed'] = StringIO('chr1\t5\t20\texon1\t.\t-')
        slices = get_slices(params)
        self.assertEqual(expected_bed, slices.head(as_string=True))
        self.assertEqual(expected_fasta, slices.print_sequence())


if __name__ == '__main__':
    unittest.main()
