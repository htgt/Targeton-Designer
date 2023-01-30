# Targeton Designer

Standalone targeton designer tool.

[[_TOC_]]

## Installation

Dependencies:

Build-essential, BedTools and Python (3.8), Python-venv (3.8)
Change ```python``` command to point to Python (3.8), ubuntu expects python3 to be a specific version for compatibility.
```sh
sudo apt-get update \
&& sudo apt-get -y install build-essential bedtools python3.8-dev python3.8-venv \
&& sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.8 2  \
&& sudo update-alternatives --config python
```

### Python3

Check Python3 (base) and Python (updated) version
```sh
python3 --version
python --version
```

### Clone the repo
Pull down the Targeton Designer repo and cd into it.
Recursively pull any submodules.
```sh
git clone --recurse-submodule https://gitlab.internal.sanger.ac.uk/sci/targeton-designer.git
cd targeton-designer
```

### Python Virtual Environment

Requirements:
 - Python3.8+
 - Python-venv

Setting up Virtual Env:
```sh
python -m venv venv

source venv/bin/activate

pip install -U pip wheel setuptools 
pip install -r requirements.txt
pip install -r sge-primer-scoring/requirements.txt

deactivate
```

Run the tests:
```sh
source venv/bin/activate

python -m unittest
cd sge-primer-scoring
python -m unittest

deactivate
```

### Exonerate iPCRess

Nathan Weeks has placed Exonerate onto GitHub along with maintenance tweaks. The repo can be found here:
https://github.com/nathanweeks/exonerate

Install glib and autoconf
```sh
sudo apt-get install libglib2.0-dev
sudo apt install autoconf
```

Installing exonerate on your VM.
```sh
git clone https://github.com/nathanweeks/exonerate.git
cd exonerate
autoreconf -vfi
./configure [YOUR_CONFIGURE_OPTIONS]
make
make check
sudo make install
```

If you happen to get the error: no module found 'apt_pkg', create a symbolic link to your apt_pkg.so

```sh
cd /usr/lib/python3/dist-packages
sudo ln -s apt_pkg.cpython-{version-number}-x86_64-linux-gnu.so apt_pkg.so
```

### Docker images

Upcoming feature in later releases

## Usage

### Command Line

Make designer.sh executable
```sh
chmod +x ./designer.sh
```

Check Designer Version:
```sh
./designer.sh version
```

#### Designer Workflow

Running full Designer Workflow:
```sh
./designer.sh design [-h] [--bed INPUT_BED] [--fasta REF_FASTA]
```

Example Command
```sh
./designer.sh design --bed example.bed --fasta example_genomic_ref.fasta
```

#### Slicer Tool

Running Slicer tool:
```sh
./designer.sh slicer [-h] [-f5 FLANK_5] [-f3 FLANK_3] [-l LENGTH] [-o OFFSET] [-d DIR] [--bed INPUT_BED] [--fasta REF_FASTA]
```

Example command:
```sh
./designer.sh slicer --bed example.bed --fasta example_genomic_ref.fa -d example_dir
```

#### Primer3 Runner

To set up user defined configuration for primer3, you can create file named ``primer3.config.json`` inside of the ``/config`` folder. 
You can use ``/config/primer3_example.config.json`` file as a template. 

If there is no config file for primer3, it will run with the default configuration.

Running Primer3:
```sh
./designer.sh primer [--fasta INPUT_FASTA] [--dir OUTPUT_FOLDER] 
```
The input fasta and BED files are intended to be sourced from the slicer tool. Examples of how these files are constructed can be found below.

Example command:
```sh
./designer.sh primer --fasta slices.fa --dir p3_output
```

#### Exonerate iPCRess

Running Exonerate iPCRess:
```sh
./designer.sh ipcress [--dir DIR] [--fasta REF_FASTA] [--primers IPCRESS_INPUT] [--p3_csv PRIMER3_OUTPUT_CSV] 
```
Supply either a standard iPCRess input file or point P3 CSV to the output csv of the Primer3 runner.

Example command:
```sh
./designer.sh ipcress --dir example_dir --fasta example_genomic_ref.fa --primers example_ipcress_input.txt
or
./designer.sh ipcress --dir example_dir --fasta example_genomic_ref.fa --p3_csv example_p3_output.csv
```

#### Primer scoring

Running primer scoring:
```sh
./designer.sh scoring [--ipcress_file IPCRESS_FILE] [--scoring_mismatch SCORING_MISMATCH] [--output_tsv OUTPUT_TSV] [--targeton_csv TARGETON_CSV] 
```

Example command:
```sh
./designer.sh scoring --ipcress_file example_ipcress_file.txt --scoring_mismatch 4 --output_tsv example_output.tsv
```
Example command with targeton csv:
```sh
./designer.sh scoring --ipcress_file example_ipcress_file.txt --scoring_mismatch 4 --output_tsv example_targeton_output.tsv --targeton_csv example_targetons.csv
```

For more information and example files see the [Primer Scoring repo](https://gitlab.internal.sanger.ac.uk/sci/sge-primer-scoring).

### Docker

Upcoming feature in later releases

## File formats
### Genomic Reference file
A Fasta file of latest GRCh38 genome. This is used for gathering the slice sequences and retrieving primer information. 
Either supply a local genome reference file or download one from EnsEMBL and point to it with the relevant parameters:
http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/

We've included a bash script to download the reference genome:
```sh
./download_reference_genome.sh
```

Otherwise you can download it manually here:
```sh
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
```

### Slicer Input BED File
A BED file containing the regions you wish to slice across. 

The chromosome column data must match your reference fasta file IDs. If youre reference had >chr1 then you must call chromosome 1 'chr1' in this column and vice-versa.

Note: BED effectively are applied tsv files so use tabs to separate the values. Headers are optional in BED file and can be a cause of issues if they aren't perfect. Strand is required for the slicer to ensure sequences are output in the correct orientation. Score isn't used but the field must be present for the file format to be read correctly.

| chrom | chromStart | chromEnd | name | score | strand |
| ----- | ---------- | -------- | ---- | ----- | ------ |
| 1 | 42931046 | 42931206 | ENSE00003571441_HG6	| 0	| - |
| 1	| 42929593 | 42929780 | ENSE00000769557_HG8 | 0 | - |

Raw file
```
1	42931046	42931206	ENSE00003571441_HG6	0	-
1	42929593	42929780	ENSE00000769557_HG8	0	-
```

More information can be found here: https://en.wikipedia.org/wiki/BED_(file_format)

### Slicer BED output
BED file output with row for each slice. This file will also be used for running VaLiAnt.

| chrom | chromStart | chromEnd | name | score | strand |
| ----- | ---------- | -------- | ---- | ----- | ------ |
| 1 | 42930996 | 42931206 | ENSE00003571441_HG6_1 | 0 | - | 
| 1 | 42931001 | 42931211 | ENSE00003571441_HG6_2 | 0 | - | 
| 1 | 42931006 | 42931216 | ENSE00003571441_HG6_3 | 0 | - |
| 1 | 42931011 | 42931221 | ENSE00003571441_HG6_4 | 0 | - |
| 1 | 42931016 | 42931226 | ENSE00003571441_HG6_5 | 0 | - |

Raw file
```
1	42930996	42931206	ENSE00003571441_HG6_1	0	-
1	42931001	42931211	ENSE00003571441_HG6_2	0	-
1	42931006	42931216	ENSE00003571441_HG6_3	0	-
1	42931011	42931221	ENSE00003571441_HG6_4	0	-
1	42931016	42931226	ENSE00003571441_HG6_5	0	-
```

### Primer3 Fasta Input File (Slicer Fasta output)
Contains the slice sequences, with their IDs including an increment, coordinates and strand in the header

```
>ENSE00003571441_HG6_1::1:42930996-42931206(-)
GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCC
TTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGCCCTCCT
>ENSE00003571441_HG6_2::1:42931001-42931211(-)
TGCAGGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTG
GCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGCC
>ENSE00003571441_HG6_3::1:42931006-42931216(-)
CCCCCTGCAGGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCAT
GATTGGCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCC
>ENSE00003571441_HG6_4::1:42931011-42931221(-)
CATCTCCCCCTGCAGGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGG
GGCATGATTGGCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAG
>ENSE00003571441_HG6_5::1:42931016-42931226(-)
GGCTGCATCTCCCCCTGCAGGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTG
TTGGGGGCATGATTGGCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCT
```

### Primer3 Output BED file
Genomic locations of the primers and their names. Names are incremented from 0 and given F and R depending on whether they're 5' or 3'

| chrom | chromStart | chromEnd | name | score | strand |
| ----- | ---------- | -------- | ---- | ----- | ------ |
| 1 | 42931021 | 42931039 | ENSE00003571441_HG6_6_LibAmpR_0 | 0 | - |
| 1 | 42931210 | 42931230 | ENSE00003571441_HG6_6_LibAmpF_0 | 0 | - |
| 1 | 42931021 | 42931039 | ENSE00003571441_HG6_6_LibAmpR_1 | 0 | - | 
| 1 | 42931211 | 42931230 | ENSE00003571441_HG6_6_LibAmpF_1 | 0 | - |

Raw file
```
1	42931021	42931039	ENSE00003571441_HG6_6_LibAmpR_0	0	-
1	42931210	42931230	ENSE00003571441_HG6_6_LibAmpF_0	0	-
1	42931021	42931039	ENSE00003571441_HG6_6_LibAmpR_1	0	-
1	42931211	42931230	ENSE00003571441_HG6_6_LibAmpF_1	0	-
```

### Primer3 Output CSV file (iPCRess Primer3 CSV Input file)
Contains all of the extra information from Primer3 for the individual primers

| primer | sequence | chr | primer_start | primer_end | tm | gc_percent | penalty | self_any_th | self_end_th | hairpin_th | end_stability |
| ------ | -------- | --- | ------------ | ---------- | -- | ---------- | ------- | ----------- | ----------- | ---------- | ------------- |
| ENSE00003571441_HG6_6_LibAmpR_0 | ACCCAGGCTGCATCTCCC | 1 | 42931021 | 42931039 | 61.41508744063151 | 66.66666666666667 | 3.4150874406315097 | 9.564684449038168 | 0.0 | 0.0 | 4.3 |
| ENSE00003571441_HG6_6_LibAmpF_0 | AGTGCCAGGACCTCTCCTAC | 1 | 42931210 | 42931230 | 60.32483047348552 | 60.0 | 0.32483047348551963 | 0.0 | 0.0 | 46.300612411542886 | 3.18 |
| ENSE00003571441_HG6_6_LibAmpR_1 | ACCCAGGCTGCATCTCCC | 1 | 42931021 | 42931039 | 61.41508744063151 | 66.66666666666667 | 3.4150874406315097 | 9.564684449038168 | 0.0 | 0.0 | 4.3 | 
| ENSE00003571441_HG6_6_LibAmpF_1 | AGTGCCAGGACCTCTCCTA | 1 | 42931211 | 42931230 | 58.90293358584404 | 57.89473684210526 | 2.097066414155961 | 0.0 | 0.0 | 46.300612411542886 | 2.94 | 

Raw File
```
primer,sequence,chr,primer_start,primer_end,tm,gc_percent,penalty,self_any_th,self_end_th,hairpin_th,end_stability
ENSE00003571441_HG6_6_LibAmpR_0,ACCCAGGCTGCATCTCCC,1,42931021,42931039,61.41508744063151,66.66666666666667,3.4150874406315097,9.564684449038168,0.0,0.0,4.3
ENSE00003571441_HG6_6_LibAmpF_0,AGTGCCAGGACCTCTCCTAC,1,42931210,42931230,60.32483047348552,60.0,0.32483047348551963,0.0,0.0,46.300612411542886,3.18
ENSE00003571441_HG6_6_LibAmpR_1,ACCCAGGCTGCATCTCCC,1,42931021,42931039,61.41508744063151,66.66666666666667,3.4150874406315097,9.564684449038168,0.0,0.0,4.3
ENSE00003571441_HG6_6_LibAmpF_1,AGTGCCAGGACCTCTCCTA,1,42931211,42931230,58.90293358584404,57.89473684210526,2.097066414155961,0.0,0.0,46.300612411542886,2.94
```

### iPCRess Standard Input File
Standard file format input for iPCRess. If you wish to use iPCRess as a standalone, use this file format.
Otherwise point iPCRess at the Primer3 output CSV.
Space separated text file containing name, left & right primer, minimum & maximum amplicon lengths.

| name | left | right | min_amplicon_len | max_amplicon_len |
| ---- | ---- | ----- | ---------------- | ---------------- |
| ENSE00003571441_HG6_6_LibAmp_0 | AGTGCCAGGACCTCTCCTAC | GGGAGATGCAGCCTGGGT | 200 | 300 |

Raw file
```
ENSE00003571441_HG6_6_LibAmp_0 AGTGCCAGGACCTCTCCTAC GGGAGATGCAGCCTGGGT 200 300
```

### iPCRess Output example
Space separated text file. Sequence_id contains the chromosome and description can be 'forward', 'revcomp', 'single_A' or 'single_B'.
| sequence_id | experiment_id | product_length | primer_5 | pos_5 | mismatch_5 | primer_3 | pos_3 | mismatch_3 | description |
| ----------- | ------------- | -------------- | -------- | ----- | ---------- | -------- | ----- | ---------- | ----------- |
| 19:filter(unmasked) | ID0001 | 259 | A | 44907726 | 0 | B | 44907967 | 0 | forward |

Raw file
```
ipcress: 19:filter(unmasked) ID0001 259 A 44907726 0 B 44907967 0 forward
```
