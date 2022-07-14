import argparse

class ParsedInputArguments:
    def __init__(self) -> None:
        self.arguments = []
        self.command = ''

        self.parse_arguments()

    def parse_arguments(self):
        parser = argparse.ArgumentParser(
            description='Targeton Designer CLI')

        parser.add_argument('command',
                            help='Command to run in Designer CLI, available commands: version, slicer, primer3, ipcress, design',
                            type=str, default='design')

        parser = add_input_args(parser)

        self.set_args(vars(parser.parse_args()))

    def get_command(self) -> str:
        return self.command

    def set_args(self, values) -> None:
        self.arguments = values
        self.command = self.arguments['command']

    def get_args(self) -> str:
        return self.arguments

def positive_int(arg: str) -> int:
    if int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0')
    return int(arg)

def len_positive_int(arg: str) -> int:
    if 10000 < int(arg) or int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0 and below 10000')
    return int(arg)

def add_input_args(parser):
    parser.add_argument('--bed',
                        help='BED file containing coords, strand and IDs')
    parser.add_argument('--fasta',
                        help='FASTA file containing seqs and IDs')
    parser.add_argument('-f5', '--flank_5',
                        help='how far to extend region at 5\' end (default 50nt)',
                        type=int, default=50)
    parser.add_argument('-1b',
                        help='Declare if input BED/TSV is one-based',
                        action='store_true')
    parser.add_argument('-f3', '--flank_3',
                        help='how far to extend region at 3\' end (default 50nt)',
                        type=int, default=50)
    parser.add_argument('-l', '--length',
                        help='length of each slice (default 210nt)',
                        type=len_positive_int, default=210)
    parser.add_argument('-o', '--offset',
                        help='offset between each slice (default 5nt)',
                        type=positive_int, default=5)
    parser.add_argument('-d', '--dir',
                        help='output directory name to be timestamped'
                             ' (default \'td_output\')',
                        type=str, default='td_output')
    parser.add_argument('--min',
                        help='Minimum amplicon length',
                        default='200')
    parser.add_argument('--max',
                        help='Maximum amplicon length',
                        default='300')
    parser.add_argument('--mismatch',
                        help='Number of mismatches to check against',
                        default='5')
    parser.add_argument('--primers',
                        help='Optional: Supply a preformatted txt file.\nIf left blank, the runner will take the primer3 output csv. Either primers or p3_csv must be supplied.')
    parser.add_argument('--p3_csv',
                        help='Optional: Point at specific Primer3 output CSV file. Either primers or p3_csv must be supplied.')
    parser.add_argument('--pretty',
                        help='Specify to include graphs in the iPCRess output. Default: false',
                        action='store_true')

    return parser

    