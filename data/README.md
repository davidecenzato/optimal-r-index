# Dataset Download

The raw sequencing data for the datasets used in the experiments are listed in `data.txt` and can be downloaded using the NCBI SRA accession numbers provided for each one.

## File format

The file `data.txt` lists the 10 datasets. For each dataset, the file provides:

1. The **dataset name**
2. A list of **NCBI SRA accession numbers**, separated by commas (e.g. `ERR732065,ERR732066`).

Datasets with more than one accession number are obtained by concatenating the corresponding files.
Some datasets list only one accession number, optionally followed by `_1`; this indicates that only the first file of a paired-end collection should be used.

## Downloading the data

You can use the `fastq-dump` command from the [SRA Toolkit](https://github.com/ncbi/sra-tools).

For each accession number, run:

```bash
fastq-dump --fasta --split-3 <ACCESSION>
```

- `--fasta` outputs the sequences in FASTA format.
- `--split-3` splits paired-end reads into separate `_1` and `_2` files, which is necessary for correctly handling paired-end sequencing runs.

For datasets with multiple accession numbers, download each accession separately, then concatenate the resulting files into a single FASTA file.
For datasets whose accession number includes `_1`, use only the corresponding `_1` file.
