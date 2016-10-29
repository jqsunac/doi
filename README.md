# PROTOCOL

Date: 2016-10-26


## 1 Preparation

### 1.1 Data preparation
Working directories and Ig-Seq data should be prepared before analysis.
The working directories and analysis scripts can be downloaded with `git` command,
and the Ig-Seq data can be downloaded from DDBJ with `wget` command.

```bash
cd ~/Desktop

## download anlaysis scripts, database, etc.
git clone git@github.com:jqsunac/fugu.git

## download Ig-Seq data.
wget http://path_to_DRR004021_1.fastq.bz2 -p ./data
wget http://path_to_DRR004021_2.fastq.bz2 -p ./data
wget http://path_to_DRR004022_1.fastq.bz2 -p ./data
wget http://path_to_DRR004022_2.fastq.bz2 -p ./data
wget http://path_to_DRR004023_1.fastq.bz2 -p ./data
wget http://path_to_DRR004023_2.fastq.bz2 -p ./data

## Decompress Ig-Seq data to FASTQ format.
bzip2 -dc ./data/DRR004021_1.fastq.bz2 > ./work/fugu1_1.fq
bzip2 -dc ./data/DRR004021_2.fastq.bz2 > ./work/fugu1_2.fq
bzip2 -dc ./data/DRR004022_1.fastq.bz2 > ./work/fugu2_1.fq
bzip2 -dc ./data/DRR004022_2.fastq.bz2 > ./work/fugu2_2.fq
bzip2 -dc ./data/DRR004023_1.fastq.bz2 > ./work/fugu3_1.fq
bzip2 -dc ./data/DRR004023_2.fastq.bz2 > ./work/fugu3_2.fq
```

| Fugu   | Number of reads |
|--------|----------------:|
| fugu 1 |       5,881,846 |
| fugu 2 |       6,199,423 |
| fugu 3 |       5,936,417 |


### 1.2 Software preparation

Then, install several software to analysis Ig-Seq data.

```bash
brew install blast
brew install fastqc
pip install cutadapt
pip install numpy
pip install matplotlib
pip install pandas
pip install biopython
pip install pydair
```

In addition, [PEAR](http://sco.h-its.org/exelixis/web/software/pear/),
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), and
[CD-HIT](http://weizhongli-lab.org/cd-hit/) also should
be installed. The following table shows that software name and version that used in this study.

| Software       | Version  | Reference                                                 |
|----------------|---------:|-----------------------------------------------------------|
| Trimmomatic    |     0.35 | http://www.usadellab.org/cms/?page=trimmomatic            |
| PEAR           |    0.9.6 | http://sco.h-its.org/exelixis/web/software/pear/          |
| CD-HIT         |   
| cutadatp       |     1.11 | http://cutadapt.readthedocs.io/en/stable/index.html       |
| PyDAIR         |    0.1.8 | https://github.com/bioinfoteam/PyDAIR                     |
| NCBI BLAST+    |    2.4.0 | https://www.ncbi.nlm.nih.gov/books/NBK279671/             |
| FastQC         |   0.11.5 | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| R              |    3.2.4 | https://cran.r-project.org/                               |



## 2 Read Classification

Each read in FASTQ file contains two types of primers that are V primer and C primer.
In this step, we sort reads according the V and C primers using [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html).

Identify the primers that used in each read.

```bash
pwd
## ~/Desktop/fugu/work

# from fugu1 to fugu4
for (( i = 1; i < 4; ++i ))
do
    # from *_1 (forward) to *_2 (reversed) of FASTQ files
    for (( j = 1; j < 3; ++j ))
    do
        cutadapt -g  nFVH1=CTGACCCAGTCTGAACCAGT   \
                 -g  nFVH2=TGAACAGTTGACACAGCCAGC  \
                 -g  nFVH3=GCCTGAAGTAAAAAGACCTGGA \
                 -g nVhCm1=CGTTCATGGTTGGAGGGTAC   \
                 -g nVhCt1=TCTGGGAAGAAGTCGAGAGC   \
                 --info-file fugu${i}_${j}.primers.info.txt \
                 --untrimmed-output /dev/null -o /dev/null  \
                 -O 10 -e 0.2 fugu${i}_${j}.fq >> log.cutadapt.txt
    done
done
```

Sort read according to V and C primers.

```bash
pwd
## ~/Desktop/fugu/work
for (( i = 1; i < 4; ++i ))
do
    python ../bin/read_classify.py --fq1 fugu${i}_1.fq --fq2 fugu${i}_2.fq \
                                   --log1 fugu${i}_1.primers.info.txt      \
                                   --log2 fugu${i}_2.primers.info.txt
done
```

The following table shows the number of reads after sorting.
The read that only contains one primer or no primer was deleted.

| Fugu   | group | V1        | V2        | V3    | total     |
|--------|-------|----------:|----------:|------:|----------:|
| fugu 1 | Cm    | 2,704,868 |   354,384 | 4,448 | 3,063,700 |
| fugu 1 | Ct    |   646,325 | 1,767,066 | 2,750 | 2,416,141 |
| fugu 2 | Cm    | 3,664,210 |   479,558 | 6,933 | 4,150,701 |
| fugu 2 | Ct    |   297,464 | 1,422,023 | 1,945 | 1,721,432 |
| fugu 3 | Cm    | 3,186,080 |   344,051 | 4,379 | 3,534,510 |
| fugu 3 | Ct    |   227,500 | 1,873,206 | 1,114 | 2,101,820 |



## 3 Quality Control

Removed the low quality reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

<!--

Quality check by FastQC.

```bash
pwd
## ~/Desktop/fugu/work

mkdir fqreports

for (( i = 1; i < 4; ++i ))
do
    for (( j = 1; j < 2; ++j))
    do
        fastqc -o fqreports -q fugu${i}_${j}.v1.cm.fq
        fastqc -o fqreports -q fugu${i}_${j}.v2.cm.fq
        fastqc -o fqreports -q fugu${i}_${j}.v3.cm.fq
        fastqc -o fqreports -q fugu${i}_${j}.v1.ct.fq
        fastqc -o fqreports -q fugu${i}_${j}.v2.ct.fq
        fastqc -o fqreports -q fugu${i}_${j}.v3.ct.fq
    done
done
```
-->


```bash
cgene=("cm" "ct")
vgene=("v1" "v2" "v3")

for (( i = 1; i < 4; ++i ))
do
    for (( c = 0; c < ${#cgene[@]}; ++c ))
    do
        for (( v = 0; v < ${#vgene[@]}; ++v ))
        do
            java -jar ../bin/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 -threads 2 \
                 -trimlog log.trimmomatic.fugu${i}.${vgene[${v}]}.${cgene[${c}]}.txt      \
                 fugu${i}_1.${vgene[${v}]}.${cgene[${c}]}.fq \
                 fugu${i}_2.${vgene[${v}]}.${cgene[${c}]}.fq \
                 fugu${i}_1.${vgene[${v}]}.${cgene[${c}]}.qc.fq \
                 fugu${i}_1.${vgene[${v}]}.${cgene[${c}]}.qc_unpaired.fq \
                 fugu${i}_2.${vgene[${v}]}.${cgene[${c}]}.qc.fq \
                 fugu${i}_2.${vgene[${v}]}.${cgene[${c}]}.qc_unpaired.fq \
                 LEADING:10 TRAILING:10 MINLEN:30 >> log.qc.txt 2>&1
        done
    done
done
```

The following table shows the number of reads fugu after performing Trimmomatic.

| Fugu   | group | V1        | V2        | V3    | total     |
|--------|-------|----------:|----------:|------:|----------:|
| fugu 1 | Cm    | 2,704,868 |   354,384 | 4,448 | 3,063,700 |
| fugu 1 | Ct    |   646,325 | 1,767,066 | 2,750 | 2,416,141 |
| fugu 2 | Cm    | 3,664,210 |   479,558 | 6,933 | 4,150,701 |
| fugu 2 | Ct    |   297,464 | 1,422,023 | 1,945 | 1,721,432 |
| fugu 3 | Cm    | 3,186,080 |   344,051 | 4,379 | 3,534,510 |
| fugu 3 | Ct    |   227,500 | 1,873,206 | 1,114 | 2,101,820 |



## 4 Merge Paired-end Reads

Merge the paired-end reads into a consensus read using 
[PEAR](http://sco.h-its.org/exelixis/web/software/pear/doc.html).

The minimum overlap size was specified to 10 bp (`-v 10`),
the p-value for the statistical test used in PEAR was specified to 0.05 (`-p 0.05`),
and the minimum length of merged read was specified to 300 (`-n 300`).

```bash
pwd
## ~/Desktop/fugu/work
cp ../bin/pear-0.9.6-bin-64/pear-0.9.6-bin-64 ../bin/pear

cgene=("cm" "ct")
vgene=("v1" "v2" "v3")

for (( i = 1; i < 4; ++i ))
do
    for (( c = 0; c < ${#cgene[@]}; ++c ))
    do
        for (( v = 0; v < ${#vgene[@]}; ++v ))
        do
            ../bin/pear -f fugu${i}_1.${vgene[${v}]}.${cgene[${c}]}.qc.fq \
                        -r fugu${i}_2.${vgene[${v}]}.${cgene[${c}]}.qc.fq \
                        -o fugu${i}.${vgene[${v}]}.${cgene[${c}]}.pear.fq \
                        -j 2 -v 10 -n 300 -p 0.05 >> log.pear.txt
        done
    done
done


## change merged FASTQ file to FASTA file.
for (( i = 1; i < 4; ++i ))
do
    awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v1.cm.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v1.cm.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v2.cm.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v2.cm.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v3.cm.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v3.cm.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v1.ct.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v1.ct.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v2.ct.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v2.ct.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v3.ct.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v3.ct.fa
done
```


The following table shows the number of reads after performing PEAR.

| Fugu   | group | V1        | V2        | V3    | total     |
|--------|-------|----------:|----------:|------:|----------:|
| fugu 1 | Cm    | 1,464,748 |   199,911 | 3,243 | 1,667,902 |
| fugu 1 | Ct    |   541,349 | 1,584,343 | 2,405 | 2,128,097 |
| fugu 2 | Cm    | 1,882,427 |   275,222 | 4,777 | 2,162,426 |
| fugu 2 | Ct    |   243,096 | 1,249,185 | 1,696 | 1,493,977 |
| fugu 3 | Cm    | 1,679,695 |   196,113 | 3,173 | 1,878,981 |
| fugu 3 | Ct    |   184,559 | 1,693,062 |   962 | 1,878,583 |




## 5 Ig-Seq Data Analysis

### 5.1 Identification of VDJ genes

Identification of V, D, and J genes using [PyDAIR](https://pypi.python.org/pypi/PyDAIR).
BLAST with the following parameters is used to align IGH sequence to V, D, J gene databases
for identification of V, D and J genes in IGH sequence.

| Gene | match score | mismatch score | gap open penalty | gap extend penalty| word size | cutoff e-value |
|---|---:|---:|---:|---:|----:|-------:|
| V |  3 | -3 |  6 |  6 |  10 |  1e-10 |
| D |  1 | -1 |  0 |  2 |   4 |  1e-02 |
| J |  3 | -3 |  6 |  6 |   7 |  1e-05 |


```bash
pwd
## ~/Desktop/fugu/work

cgene=("m" "t")         # cm, ct
vgene=("1" "2" "3")     # v1, v2, v3

# Identify V, D, and J genes

for (( i = 1; i < 4; ++i ))
do
    for (( c = 0; c < ${#cgene[@]}; ++c ))
    do
        for (( v = 0; v < ${#vgene[@]}; ++v ))
        do
            pydair parse -q fugu${i}.v${vgene[${v}]}.c${cgene[${c}]}.fa \
                         -v ../db/V${vgene[${v}]}.fa \
                         -d ../db/D${cgene[${c}]}.fa \
                         -j ../db/J${cgene[${c}]}.fa \
                         --v-blastdb ../db/v${vgene[${v}]}db \
                         --d-blastdb ../db/d${cgene[${c}]}db \
                         --j-blastdb ../db/j${cgene[${c}]}db \
                         -o fugu${i}.v${vgene[${v}]}.c${cgene[${c}]}
        done
    done
done



## merge V1, V2 and V3 groups into one group.

cgene=("m" "t")         # cm, ct
for (( i = 1; i < 4; ++i ))
do
    for (( c = 0; c < ${#cgene[@]}; ++c ))
    do
        cat fugu${i}.v1.c${cgene[${c}]}.vdj.pydair >  fugu${i}.c${cgene[${c}]}.pydair
        cat fugu${i}.v2.c${cgene[${c}]}.vdj.pydair >> fugu${i}.c${cgene[${c}]}.pydair
        cat fugu${i}.v3.c${cgene[${c}]}.vdj.pydair >> fugu${i}.c${cgene[${c}]}.pydair
    done
done
```

### 5.2 Calculation of frequencies of VDJ usage

```bash
cgene=("m" "t")
for ((c = 0; c < ${#cgene[@]}; ++c))
do
    pydair stats -i fugu1.c${cgene[${c}]}.pydair fugu2.c${cgene[${c}]}.pydair fugu3.c${cgene[${c}]}.pydair \
                 -n Fugu1 Fugu2 Fugu3 \
                 -o fugustats_c${cgene[${c}]} \
                 --contain_ambiguous_D
done
```



### 5.3 CDR3 sequence

Create CDR3 amino acid sequence.

```bash
cgene=("m" "t")         # cm, ct
for (( i = 1; i < 4; ++i ))
do
    for (( c = 0; c < ${#cgene[@]}; ++c ))
    do
        python ../bin/create_cdr3_fasta.py -i fugu${i}.c${cgene[${c}]}.pydair  \
                                           -o fugu${i}.c${cgene[${c}]}.cdr3aa.fa \
                                           -f fugu${i}
    done
done
```


## 6 Statsitics Analysis

### 6.1 Data files

- Fugu 1
    - [Fugu 1 Ct V](work/fugustats_ct.sample.Fugu1.v.freq.tsv)
    - [Fugu 1 Ct D](work/fugustats_ct.sample.Fugu1.d.freq.tsv)
    - [Fugu 1 Ct J](work/fugustats_ct.sample.Fugu1.j.freq.tsv)
    - [Fugu 1 Ct VDJ](work/fugustats_ct.sample.Fugu1.vdj.freq.tsv)
    - [Fugu 1 Ct CDR3 AA length](work/fugustats_ct.sample.Fugu1.cdr3_prot_length.freq.tsv)
    - [Fugu 1 Ct CDR3 AA sequence](work/fugu1.ct.cdr3aa.fa.gz)
    - [Fugu 1 Cm V](work/fugustats_cm.sample.Fugu1.v.freq.tsv)
    - [Fugu 1 Cm D](work/fugustats_cm.sample.Fugu1.d.freq.tsv)
    - [Fugu 1 Cm J](work/fugustats_cm.sample.Fugu1.j.freq.tsv)
    - [Fugu 1 Cm VDJ](work/fugustats_cm.sample.Fugu1.vdj.freq.tsv)
    - [Fugu 1 Cm CDR3 AA length](work/fugustats_cm.sample.Fugu1.cdr3_prot_length.freq.tsv)
    - [Fugu 1 Cm CDR3 AA sequence](work/fugu1.cm.cdr3aa.fa.gz)
- Fugu 2
    - [Fugu 2 Ct V](work/fugustats_ct.sample.Fugu2.v.freq.tsv)
    - [Fugu 2 Ct D](work/fugustats_ct.sample.Fugu2.d.freq.tsv)
    - [Fugu 2 Ct J](work/fugustats_ct.sample.Fugu2.j.freq.tsv)
    - [Fugu 2 Ct VDJ](work/fugustats_ct.sample.Fugu2.vdj.freq.tsv)
    - [Fugu 2 Ct CDR3 AA length](work/fugustats_ct.sample.Fugu2.cdr3_prot_length.freq.tsv)
    - [Fugu 2 Ct CDR3 AA sequence](work/fugu2.ct.cdr3aa.fa.gz)
    - [Fugu 2 Cm V](work/fugustats_cm.sample.Fugu2.v.freq.tsv)
    - [Fugu 2 Cm D](work/fugustats_cm.sample.Fugu2.d.freq.tsv)
    - [Fugu 2 Cm J](work/fugustats_cm.sample.Fugu2.j.freq.tsv)
    - [Fugu 2 Cm VDJ](work/fugustats_cm.sample.Fugu2.vdj.freq.tsv)
    - [Fugu 2 Cm CDR3 AA length](work/fugustats_cm.sample.Fugu2.cdr3_prot_length.freq.tsv)
    - [Fugu 2 Cm CDR3 AA sequence](work/fugu2.cm.cdr3aa.fa.gz)
- Fugu 3
    - [Fugu 3 Ct V](work/fugustats_ct.sample.Fugu3.v.freq.tsv)
    - [Fugu 3 Ct D](work/fugustats_ct.sample.Fugu3.d.freq.tsv)
    - [Fugu 3 Ct J](work/fugustats_ct.sample.Fugu3.j.freq.tsv)
    - [Fugu 3 Ct VDJ](work/fugustats_ct.sample.Fugu3.vdj.freq.tsv)
    - [Fugu 3 Ct CDR3 AA length](work/fugustats_ct.sample.Fugu3.cdr3_prot_length.freq.tsv)
    - [Fugu 3 Ct CDR3 AA sequence](work/fugu3.ct.cdr3aa.fa.gz)
    - [Fugu 3 Cm V](work/fugustats_cm.sample.Fugu3.v.freq.tsv)
    - [Fugu 3 Cm D](work/fugustats_cm.sample.Fugu3.d.freq.tsv)
    - [Fugu 3 Cm J](work/fugustats_cm.sample.Fugu3.j.freq.tsv)
    - [Fugu 3 Cm VDJ](work/fugustats_cm.sample.Fugu3.vdj.freq.tsv)
    - [Fugu 3 Cm CDR3 AA length](work/fugustats_cm.sample.Fugu3.cdr3_prot_length.freq.tsv)
    - [Fugu 3 Cm CDR3 AA sequence](work/fugu3.cm.cdr3aa.fa.gz)




### 6.2 Frequency of V, D, J Genes Usage

Plot the frequency of V, D, and J gene usages using R script for each group.

```
pwd
## ~/Desktop/fugu/work

Rscript --vanilla --slave ../bin/plot_freq.R \
          'fugustats_cm.sample.Fugu1,fugustats_cm.sample.Fugu2,fugustats_cm.sample.Fugu3' \
          'Fugu 1,Fugu 2,Fugu 3' TRUE Fig_Cm_VDJ_freq_hist.png png 800 400 100
Rscript --vanilla --slave ../bin/plot_freq.R \
          'fugustats_ct.sample.Fugu1,fugustats_ct.sample.Fugu2,fugustats_ct.sample.Fugu3' \
          'Fugu 1,Fugu 2,Fugu 3' TRUE Fig_Ct_VDJ_freq_hist.png png 800 400 100
```

---

**Cm genes**

![Cm gene usage frequency](./work/Fig_Cm_VDJ_freq_hist.png)

---

**Ct genes**

![Ct gene usage frequency](./work/Fig_Ct_VDJ_freq_hist.png)

---



### 6.3 VDJ Combinations

Plot 3D figures for visualizing VDJ combinations.

```bash
Rscript --vanilla --slave ../bin/plot_vdj_3d.R \
            fugustats_cm.sample.Fugu1.vdj.freq.tsv Fugu1 TRUE Fig_Cm_VDJ3D_fugu1
Rscript --vanilla --slave ../bin/plot_vdj_3d.R \
            fugustats_cm.sample.Fugu2.vdj.freq.tsv Fugu2 TRUE Fig_Cm_VDJ3D_fugu2
Rscript --vanilla --slave ../bin/plot_vdj_3d.R \
            fugustats_cm.sample.Fugu3.vdj.freq.tsv Fugu3 TRUE Fig_Cm_VDJ3D_fugu3
```

---

**Fugu 1**

![Cm VDJ combinations for Fugu 1](./work/Fig_Cm_VDJ3D_fugu1.3d.png)
![Cm VDJ combinations for Fugu 1](./work/Fig_Cm_VDJ3D_fugu1.color.png)

---

**Fugu 2**

![Cm VDJ combinations for Fugu 2](./work/Fig_Cm_VDJ3D_fugu2.3d.png)
![Cm VDJ combinations for Fugu 1](./work/Fig_Cm_VDJ3D_fugu2.color.png)

---

**Fugu 3**

![Cm VDJ combinations for Fugu 3](./work/Fig_Cm_VDJ3D_fugu3.3d.png)
![Cm VDJ combinations for Fugu 1](./work/Fig_Cm_VDJ3D_fugu3.color.png)


---



### 6.4 Rarefaction Study for VDJ Combinations.

Perform rarefaction study for VDJ combinations using R script for Cm group.

```
pwd
## ~/Desktop/fugu/work

Rscript --vanilla --slave ../bin/plot_rarefaction_vdjcombi.R \
          'fugustats_cm.sample.Fugu1.vdj.freq.tsv,fugustats_cm.sample.Fugu2.vdj.freq.tsv,fugustats_cm.sample.Fugu3.vdj.freq.tsv' \
          'Fugu 1,Fugu 2,Fugu 3' Fig_Cm_VDJ_rarefaction.png png 800 400 100
```

![VDJ Combinations Rarefaction Curve](./work/Fig_Cm_VDJ_rarefaction.png)



### 6.5 CDR3 Length Distributions

Plot the distributions of CDR3 amino acid sequences for each fugu data.



```
pwd
## ~/Desktop/fugu/work

Rscript --vanilla --slave ../bin/plot_cdr3len.R \
        'fugustats_cm.sample.Fugu1.cdr3_prot_length.freq.tsv,fugustats_cm.sample.Fugu2.cdr3_prot_length.freq.tsv,fugustats_cm.sample.Fugu3.cdr3_prot_length.freq.tsv' \
        'Fugu 1,Fugu 2,Fugu 3' TRUE Fig_Cm_CDR3Length.png png 800 300 100
Rscript --vanilla --slave ../bin/plot_cdr3len.R \
        'fugustats_ct.sample.Fugu1.cdr3_prot_length.freq.tsv,fugustats_ct.sample.Fugu2.cdr3_prot_length.freq.tsv,fugustats_ct.sample.Fugu3.cdr3_prot_length.freq.tsv' \
        'Fugu 1,Fugu 2,Fugu 3' TRUE Fig_Ct_CDR3Length.png png 800 300 100
```

---

**Cm group**

![Cm CDR3 Length](./work/Fig_Cm_CDR3Length.png)

---

**Ct group**

![Ct CDR3 Length](./work/Fig_Ct_CDR3Length.png)

---



### 6.6 CDR3 diversity study

CDR3 amino acid sequences of the three fugu are merged into one file for estimating CDR3 repertoire population sizes.
Then, CDR3 amino acid sequences were clustered by [CD-HIT](http://weizhongli-lab.org/cd-hit/). Output files from CD-HIT were adjusted by Python script, and estimation was performed by R script.


```bash
pwd
## ~/Desktop/fuug/work

cgene=("m" "t")

## merge fugu1, fugu2, fugu3 data into one file.
for ((c = 0; c < ${#cgene[@]}; ++c))
do
    cat fugu1.c${cgene[${c}]}.cdr3aa.fa > fuguall.c${cgene[${c}]}.cdr3aa.fa
    cat fugu2.c${cgene[${c}]}.cdr3aa.fa >> fuguall.c${cgene[${c}]}.cdr3aa.fa
    cat fugu3.c${cgene[${c}]}.cdr3aa.fa >> fuguall.c${cgene[${c}]}.cdr3aa.fa
done

## cluster CDR3 sequences using CD-HIT
for ((c = 0; c < ${#cgene[@]}; ++c))
do
    ../bin/cd-hit-v4.6.6-2016-0711/cd-hit \
             -i fuguall.c${cgene[${c}]}.cdr3aa.fa \
             -o cdhit_fuguall.c${cgene[${c}]}.cls.fa \
             -c 0.8 -n 4 -s 0.8 -M 2000 -l 5 -d 200
done

## create dataset for diversity study
for ((c = 0; c < ${#cgene[@]}; ++c))
do
    ../bin/make_cdr3_capturerecapture_data.py \
             -i cdhit_fuguall.c${cgene[${c}]}.cls.fa.clstr \
             -o fuguall.c${cgene[${c}]}.capturestudy.dataset.tsv \
             -p all
done

## CDR3 diversity sutdy
for ((c = 0; c < ${#cgene[@]}; ++c))
do
    Rscript ../bin/study_cdr3aaseq.R c${cgene[${c}]} fugustatscdr3
done
```

The following table shows the number of clusters that clustered by CD-HIT for each fugu.

| Fugu   | Cm        |Ct        |
|--------|----------:|---------:|
| Fugu 1 |   166,469 |   28,566 |
| Fugu 2 |   212,230 |   23,031 |
| Fugu 3 |   191,320 |   30,934 |


The figures shows the overlaps of the clusters of CDR3 amino acid sequences among the three fugu of Cm group.


![Cm CDR3 Venn](./work/fugustatscdr3.CDR3aaseq.cm.venn.png)



The figures shows the overlaps of the clusters of CDR3 amino acid sequences among the three fugu of Ct group.


![Ct CDR3 Venn](./work/fugustatscdr3.CDR3aaseq.ct.venn.png)

The following table shows the estimated size of CDR3 repertore.


|                | Cm      | Ct      |
|----------------|--------:|--------:|
| Estimated Size | 617,929 | 145,627 |


Create TSV files for showing the common sequences. The TSV file consists of 8 columns:

- [CDR3 clusters Cm](./work/fuguall.cm.cdr3.venn.tsv.gz)
- [CDR3 clusters Ct](./work/fuguall.ct.cdr3.venn.tsv.gz)

1. Cluster ID
2. Representative sequence
3. Representative sequence length
4. Number of sequences in this class of Fugu 1
5. Number of sequences in this class of Fugu 2
6. Number of sequences in this class of Fugu 3
7. Number of sequences in this class (Fugu1 + Fugu2 + Fugu3) 
8. Sequence IDs

```bash
python ../bin/make_cdr3_clusters_tsv.py -f fuguall.cm.cdr3aa.fa          \
                                        -c cdhit_fuguall.cm.cls.fa.clstr \
                                        -o fuguall.cm.cdr3.venn.tsv
python ../bin/make_cdr3_clusters_tsv.py -f fuguall.ct.cdr3aa.fa          \
                                        -c cdhit_fuguall.ct.cls.fa.clstr \
                                        -o fuguall.ct.cdr3.venn.tsv                                        
```

