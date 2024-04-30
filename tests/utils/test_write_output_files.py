import unittest
from os import path
from pathlib import Path
import csv

import pandas as pd
from pyfakefs.fake_filesystem_unittest import TestCase
from unittest.mock import patch, Mock
from freezegun import freeze_time

from primer.write_primer_output import _add_double_quotes_to_str, _round_to_three_decimals
from src.utils.write_output_files import write_scoring_output, write_targeton_csv
from src.utils import write_output_files
from primer.slice_data import SliceData


class TestWriteOutputFiles(TestCase):
    def setUp(self):
        self.setUpPyfakefs()
        self.fs.create_dir('test_dir')
        contents = (
            'region_1_1 ATCG GCTA 200 300\n'
            'region_1_2 ATCG GCTA 200 300\n'
            'region_2_1 ATCG GCTA 200 300\n'
        )
        self.fs.create_file('test_ipcress_input.txt', contents=contents)
        contents = (
            'chr1\t100\t250\tregion_1\t.\t+\n'
            'chr2\t100\t250\tregion_2\t.\t+\n'
        )
        self.fs.create_file('test.bed', contents=contents)
        self.slices = [
            SliceData('region_1', 'start', 'end', 'strand', 'chrom', 'bases'),
            SliceData('region_2', 'start', 'end', 'strand', 'chrom', 'bases'),
        ]

    @patch('builtins.print')
    def test_write_targeton_csv_success_with_timestamped_dir(self, mock_print):
        # arrange
        self.fs.create_dir('test_dir/td_310123')
        expected_dir = 'test_dir/td_310123'
        expected_file = 'test_dir/td_310123/targetons.csv'

        # act
        result = write_targeton_csv(
            'test_ipcress_input.txt', self.slices, 'test_dir/td_310123', dir_timestamped=True
        )

        # assert
        self.assertTrue(path.exists(expected_file))
        mock_print.assert_called_with('Targeton csv generated: test_dir/td_310123/targetons.csv')
        self.assertEqual(result.dir, expected_dir)
        self.assertEqual(result.csv, expected_file)

    @patch('builtins.print')
    @freeze_time('2023-01-31')
    def test_write_targeton_csv_success_with_non_timestamped_dir(self, mock_print):
        # arrange
        expected_dir = 'test_dir/td_20230131000000000000'
        expected_file = 'test_dir/td_20230131000000000000/targetons.csv'

        # act
        result = write_targeton_csv('test_ipcress_input.txt', self.slices, 'test_dir')

        # assert
        self.assertTrue(path.exists(expected_file))
        mock_print.assert_called_with(
            'Targeton csv generated: test_dir/td_20230131000000000000/targetons.csv'
        )
        self.assertEqual(result.dir, expected_dir)
        self.assertEqual(result.csv, expected_file)

    @patch('builtins.print')
    def test_write_targeton_csv_makes_csv_with_correct_contents(self, mock_print):
        # arrange
        expected = (
            'region_1_1,region_1\n'
            'region_1_2,region_1\n'
            'region_2_1,region_2\n'
        )

        # act
        write_targeton_csv('test_ipcress_input.txt', self.slices, 'test_dir', True)
        with open('test_dir/targetons.csv') as f:
            actual = f.read()

        # assert
        self.assertEqual(actual, expected)

    @patch('builtins.print')
    def test_write_scoring_output_success(self, mock_print):
        # arrange
        mock_scoring = Mock()
        expected_dir = ''
        expected_file = 'test_scoring.tsv'

        # act
        result = write_scoring_output(mock_scoring, 'test_scoring.tsv')

        # assert
        mock_scoring.save_mismatches.assert_called_with(expected_file)
        mock_print.assert_called_with('Scoring file saved: test_scoring.tsv')
        self.assertEqual(result.dir, expected_dir)
        self.assertEqual(result.tsv, expected_file)

    def test_export_to_csv(self):
        # arrange
        expected_file = 'test_scoring.tsv'
        fake_dir = 'fake_dir/test'
        self.fs.create_dir(fake_dir)
        mock_data = {'test': [1, 2, 3, 4, 5], 'test2': 'things'}
        mock_headers = ['test', 'test2']
        expected_delimiter = '\t'
        expected_file_path = Path(fake_dir) / expected_file

        # act
        result = write_output_files.export_to_csv(
            mock_data,
            fake_dir,
            expected_file,
            mock_headers,
            delimiter=expected_delimiter
        )

        sniffer = csv.Sniffer()
        with open(expected_file_path) as f:
            test_delimiter = sniffer.sniff(f.read()).delimiter
            f.seek(0)
            test_data = f.read()
        expected_read_data = f"test{expected_delimiter}test2\n[1, 2, 3, 4, 5]{expected_delimiter}things\n"
        # assert
        self.assertEqual(result, expected_file_path)
        self.assertEqual(test_delimiter, expected_delimiter)
        self.assertTrue(expected_file_path.exists())
        self.assertEqual(test_data, expected_read_data)


class TestDataframeMapperFunctions(unittest.TestCase):

    def test_add_double_quotes_df(self):
        df = pd.DataFrame({'strings': ['hello', 'world'], 'numbers': [3.141592653589793, 1.23456789]})

        # Apply _add_double_quotes_to_str to the 'strings' column of DataFrame
        result = df.applymap(_add_double_quotes_to_str)

        # Check if the function correctly adds double quotes to strings
        self.assertEqual(result['strings'][0], "\"hello\"")
        self.assertEqual(result['strings'][1], "\"world\"")

    def test_round_to_three_decimals_df(self):
        df = pd.DataFrame({'strings': ['hello', 'world'], 'numbers': [3.141592653589793, 1.23456789]})

        # Apply round_to_three_decimals to the 'numbers' column of DataFrame
        result = df.applymap(_round_to_three_decimals)

        # Check if the function correctly rounds floats to three decimal places
        self.assertAlmostEqual(result['numbers'][0], 3.142)
        self.assertAlmostEqual(result['numbers'][1], 1.235)

    def test_add_double_quotes_to_str(self):
        # Test adding double quotes to a string
        self.assertEqual(_add_double_quotes_to_str("hello"), "\"hello\"")
        self.assertEqual(_add_double_quotes_to_str("world"), "\"world\"")

        # Test when input is not a string
        self.assertEqual(_add_double_quotes_to_str(123), 123)

    def test_round_to_three_decimals(self):
        # Test rounding a float to three decimal places
        self.assertAlmostEqual(_round_to_three_decimals(3.141592653589793), 3.142)
        self.assertAlmostEqual(_round_to_three_decimals(1.23456789), 1.235)

        # Test when input is not a float
        self.assertEqual(_round_to_three_decimals("3.1415"), "3.1415")
        self.assertEqual(_round_to_three_decimals(123), 123)


if __name__ == '__main__':
    unittest.main()
