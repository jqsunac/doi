# PyDAIR

PyDAIR is a Python library for diversity analysis of immune repertoire.
PyDAIR provides methods for analyzing the followings features:

 - identifying the V, D, and J genes in immunoglobulin heavy chain,
   and T-cell receptor beta and delta chains.
 - identifying the CDR3 region in the sequences.
 - estimating the population sizes of the number of combinations of
   VDJ, and the number of clusters of CDR3 sequence.
 - visualizing the analysed results.

## Dependency Libraries

PyDAIR depends on the following Python libraries.

 - numpy
 - matplotlib
 - pandas
 - biopython

## Installation

Before install the PyDAIR, one need to install all depending libraries.
All dpending libraries can be installed by `pip` command.

```bash
pip install numpy
pip install matplotlib
pip install pandas
pip install biopython
``` 

Then install PyDAIR library.

```bash
pip install PyDAIR
```


# Usage Examples

## Command line usage

```sh
PyDAIR -q ./data/test.fa -o ./data/test_PyDAIR_cmd -f PyDAIR \
       -v ./data/db/v.fa -d ./data/db/d.fa -j ./data/db/j.fa \
       --v-blastdb ./data/db/v --v-match-score 3 --v-mismatch-score -3 \
       --v-gap-open-penalty 6 --v-gap-extend-penalty 6 --v-wordsize 21 --v-evalue-cutoff 1e-60 \
       --d-blastdb ./data/db/d --d-match-score 1 --d-mismatch-score -1 \
       --d-gap-open-penalty 0 --d-gap-extend-penalty 2 --d-wordsize 4 --d-evalue-cutoff 1 \
       --j-blastdb ./data/db/j --j-match-score 3 --j-mismatch-score -3 \
       --j-gap-open-penalty 6 --j-gap-extend-penalty 6 --j-wordsize 7 --j-evalue-cutoff 1e-5

```

