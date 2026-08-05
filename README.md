# optimal-r-index

This repository contains an implementation of the **optimal r-index**, i.e., the r-index of the BWT of string collections guaranteeing the minimum number of runs. Such BWT is called **optimal BWT**, and can be computed using the [optimalBWT](https://github.com/davidecenzato/optimalBWT) tool.

### Prerequesites

The optimal-r-index tool requires
* A modern Python 3 release version 3.7 or higher.
* A modern C++11 compiler such as `g++` version 4.9 or higher.
* A modern 64bits MacOs or Linux based system (MacOs not for the experimental pipeline).

### Download and Compile

```console
git clone https://github.com/davidecenzato/optimal-r-index.git
cd optimal-r-index
git submodule update --init --recursive

python3 install.py
```

### Usage

See usage options:
```console
usage: optimal-r-index.py [-h] [--algo {sais,bcr}] [--keep] [input] [output]

Tool for computing the optimal r-index of a string collection

positional arguments:
  input              input FASTA path
  output             output file path

options:
  -h, --help         show this help message and exit
  --algo {sais,bcr}  algorithm for computing the optimal BWT (sais|bcr). Default: sais
  --keep             keep temporary files
```

The current implementation takes in input string collections in FASTA format only. You can choose between two algorithms for computing the optimal BWT: **`sais`**: works in internal memory and is faster for small datasets, **`bcr`**: works in semi-external memory and is better suited for large datasets.

The `--keep` flag preserves all temporary files and is intended for debugging purposes. The main output of this software is the optimal r-index, stored in binary format as `output.ri`.

### Run on Example Data

```console
// Construct the optimal-r-index of a toy dataset
python3 optimal-r-index.py data/toy.fasta toy_index 

// Run the testing pipeline for all FASTA files in a directory
python3 exp_pipeline.py data output.csv
```

### External resources

* [optimalBWT](https://github.com/davidecenzato/optimalBWT.git)
* [BCR_LCP_GSA](https://github.com/giovannarosone/BCR_LCP_GSA.git)
* [sdsl-lite](https://github.com/simongog/sdsl-lite.git)
* [r-index](https://github.com/nicolaprezza/r-index.git)
* [Big-BWT](https://github.com/alshai/Big-BWT.git)

### Authors

* Davide Cenzato
* Veronica Guerrini
* Zsuzsanna Lipták
* Giovanna Rosone

#### Software coding and Experimental results:

* Davide Cenzato
* Veronica Guerrini

### Reference and citation 

[1] Davide Cenzato, Veronica Guerrini, Zsuzsanna Lipták, Giovanna Rosone: Computing the optimal BWT of very large string collections. DCC 2023: 71-80 ([go to the paper](https://doi.org/10.1109/DCC55655.2023.00015))

If you use this tool in an academic setting, please cite this work as follows:

**conference paper**

    @inproceedings{CenzatoGLR23,
      author       = {Davide Cenzato and
                      Veronica Guerrini and
                      Zsuzsanna Lipt{\'{a}}k and
                      Giovanna Rosone},
      title        = {Computing the optimal {BWT} of very large string collections},
      booktitle    = {In Proc. of the 33rd Data Compression Conference, {DCC} 2023},
      pages        = {71--80},
      year         = {2023},
      doi          = {10.1109/DCC55655.2023.00015}
    }

### Contacts

If you notice any bugs, please feel free to report them by opening a Git issue or by contacting us at davide_dot_cenzato_at_unive_dot_it email.
