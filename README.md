# BLIGH

BLIGH is a Python framework to analyze the diversity of immunoglobulin heavy chain.
BLIGH identifies the V, D, and J genes in the seqeunces that sequenced by High-throughput sequencer.
In addtion, BLIGH also identifies the CDR3 region.



# Usage

## Python module usage

```python
from BLIGH.seq.IgSeq import IgSeq
from BLIGH.io.BLIGHIO import *
from BLIGH.utils.BLIGHUtils import *
from BLIGH.utils.BLIGHArgs import *
from BLIGH.app.BLIGHAPP import *

# variables settings
v_gene_align_args = BLIGHBlastArgs('./data/db/v', 3, -3, 6, 6, 21, 1e-60)
d_gene_align_args = BLIGHBlastArgs('./data/db/d', 1, -1, 0, 2,  4, 1)
j_gene_align_args = BLIGHBlastArgs('./data/db/j', 3, -3, 6, 6,  7, 1e-5)
q_fasta      = './data/test.fa'
v_gene_fasta = './data/db/v.fa'
d_gene_fasta = './data/db/d.fa'
j_gene_fasta = './data/db/j.fa'
output_prefix = './data/test_bligh_moduleuse'

# BLIGH arguemnts settings
bligh_args = BLIGHArgs(q_fasta, v_gene_fasta, d_gene_fasta, j_gene_fasta,
                       output_prefix, 'bligh',
                       v_gene_align_args, d_gene_align_args, j_gene_align_args)

# main processes
bligh = BLIGHAPP(bligh_args)
bligh.blast('v')
bligh.blast('j')
bligh.parse_VJ()
bligh.write_bligh()
bligh.write_fasta('unaligned_seq')
bligh.blast('d')
bligh.parse_VDJ()
bligh.write_bligh()
```


## Command line usage

```sh
bligh -q ./data/test.fa -o ./data/test_bligh_cmd -f bligh \
      -v ./data/db/v.fa -d ./data/db/d.fa -j ./data/db/j.fa \
      --v-blastdb ./data/db/v --v-match-score 3 --v-mismatch-score -3 \
      --v-gap-open-penalty 6 --v-gap-extend-penalty 6 --v-wordsize 21 --v-evalue-cutoff 1e-60 \
      --d-blastdb ./data/db/d --d-match-score 1 --d-mismatch-score -1 \
      --d-gap-open-penalty 0 --d-gap-extend-penalty 2 --d-wordsize 4 --d-evalue-cutoff 1 \
      --j-blastdb ./data/db/j --j-match-score 3 --j-mismatch-score -3 \
      --j-gap-open-penalty 6 --j-gap-extend-penalty 6 --j-wordsize 7 --j-evalue-cutoff 1e-5

```

