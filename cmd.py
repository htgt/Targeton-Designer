import os
import sys
import csv
import argparse

from runner.primer3 import primer3_runner
from runner.exonerate import generate_ipcress_input 

parser = argparse.ArgumentParser()

parser.add_argument('--seqs', action='store', type=str, nargs='+', help='Path to sequences file')
parser.add_argument('--dir', action='store', type=str, nargs='+', help='Output folder location')

args = parser.parse_args()

file_name = args.seqs
file_csv = open(file_name[0])
csv_reader = csv.DictReader(file_csv)
rows = []
for row in csv_reader:
    rows.append(row)
file_csv.close()

run_id = rows[0]['id']

design_input = {
    'SEQUENCE_ID': run_id,
    'SEQUENCE_TEMPLATE': rows[0]['seq'],
}

dir_path = args.dir[0]
os.environ["PRIMER3_CONFIG"] = "./primer3_config.json"
primers = primer3_runner(design_input)
run_ipcress(run_id, primers, dir_path)
print(primers)
