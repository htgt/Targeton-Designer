import unittest
import sys 

from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase
from pathlib import Path
from tempfile import TemporaryDirectory

from cli import slicer_command, primer_command
from utils.arguments_parser import ParsedInputArguments


class TestSlicerIntegration(TestCase):
    def setUp(self):
        self.bed_file_path = r"./tests/integration/fixtures/bed_example.bed"
        self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"

    def test_SlicerOutput(self):
        with TemporaryDirectory() as tmpdir:
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "slicer", "--bed", self.bed_file_path, "--fasta", self.fasta_file_path, "--dir", tmpdir]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()
                slicer_result = slicer_command(args)
                path_bed = Path(slicer_result.bed)
                path_fasta = Path(slicer_result.fasta)
                # # Check if the files exist.
                self.assertTrue(path_bed.is_file())
                self.assertTrue(path_fasta.is_file())
                # # Check if the files are empty
                self.assertGreater(path_bed.stat().st_size, 0)
                self.assertGreater(path_fasta.stat().st_size, 0)


class TestPrimerIntegration(TestCase):
    def setUp(self):
        self.fasta_file_path = r"./examples/test_slice_seqs.fa"

    def test_PrimerOutput(self):
        with TemporaryDirectory() as tmpdir:
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "primer", "--fasta", self.fasta_file_path, "--dir", tmpdir]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()
                primer_result = primer_command(args["fasta"], prefix=args["dir"])
                path_bed = Path(primer_result.bed)
                path_csv = Path(primer_result.csv)
                # # Check if the files exist.
                self.assertTrue(path_bed.is_file())
                self.assertTrue(path_csv.is_file())
                # # Check if the files are empty
                self.assertGreater(path_bed.stat().st_size, 0)
                self.assertGreater(path_csv.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()