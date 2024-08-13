# optimalBWT journal paper experimental lineup

### Prerequesites

sdsl-lite library installed in /usr/local/

### Download and Compile

```console
git clone https://github.com/davidecenzato/optimalBWTJournal.git
cd optimalBWTJournal
git submodule update --init --recursive

python3 compile.py
```

### Usage

See usage options:
```console
python3 pipeline.py -h

usage: pipeline.py [-h] [--multi] [--opt] [--concat] [input_folder] [output_file]
```

input_folder: name of the folder containing the input fasta files.
output_file: name of the csv file containing the results.
--multi: compute size of the r-index using the multidollar BWT.
--opt: compute size of the r-index using the optimal BWT.
--concat: compute size of the r-index using the concatenated BWT (the original one).

Example of usage:
```console
python3 pipeline.py --multi --opt --concat data/ exp1.csv
```