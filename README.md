# PROTOCOL

Date: 2016-12-02


## 1 Preparation

### 1.1 Data preparation

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

### 1.2 Software preparation

```bash
brew install blast
brew install fastqc
brew install trimmomatic
pip install cutadapt
pip install numpy
pip install matplotlib
pip install pandas
pip install biopython
pip install pydair
```

| Software       | Version  | Reference                                                 |
|----------------|---------:|-----------------------------------------------------------|
| Trimmomatic    |     0.35 | http://www.usadellab.org/cms/?page=trimmomatic            |
| PEAR           |    0.9.6 | http://sco.h-its.org/exelixis/web/software/pear/          |
| CD-HIT         |    4.6.6 | http://weizhongli-lab.org/cd-hit/                         |
| cutadatp       |     1.11 | http://cutadapt.readthedocs.io/en/stable/index.html       |
| PyDAIR         |    0.1.8 | https://github.com/bioinfoteam/PyDAIR                     |
| NCBI BLAST+    |    2.4.0 | https://www.ncbi.nlm.nih.gov/books/NBK279671/             |
| FastQC         |   0.11.5 | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| R              |    3.2.4 | https://cran.r-project.org/                               |



## 2 Read Classification

Check the primers and split FASTQ file into (Cm, Ct) x (V1, V2, V3) groups.

```bash
pwd
## ~/Desktop/fugu/work


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


for (( i = 1; i < 4; ++i ))
do
    python ../bin/read_classify.py --fq1 fugu${i}_1.fq --fq2 fugu${i}_2.fq \
                                   --log1 fugu${i}_1.primers.info.txt      \
                                   --log2 fugu${i}_2.primers.info.txt
done
```


## 3. Quality Control



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
                 LEADING:20 TRAILING:20 MINLEN:30 >> log.qc.txt 2>&1
        done
    done
done
```



## 4 Merge Paired-end Reads

Merge the paired-end reads into a consensus read using 
[PEAR](http://sco.h-its.org/exelixis/web/software/pear/doc.html).


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



## 5 Ig-Seq Data Analysis


Identification of V, D, and J genes using [PyDAIR](https://pypi.python.org/pypi/PyDAIR).
BLAST with the following parameters is used to align IGH sequence to V, D, J gene databases
for identification of V, D and J genes in IGH sequence.


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
                         --v-blastdb ../db/v${vgene[${v}]}db --v-evalue-cutoff 1e-90 \
                         --d-blastdb ../db/d${cgene[${c}]}db \
                         --j-blastdb ../db/j${cgene[${c}]}db --j-evalue-cutoff 1e-9 \
                         -o fugu${i}.v${vgene[${v}]}.c${cgene[${c}]}
        done
    done
done

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

Summarizes the anlaysis results.


```bash
cgene=("m" "t")
for ((c = 0; c < ${#cgene[@]}; ++c))
do
    pydair stats -i fugu1.c${cgene[${c}]}.pydair fugu2.c${cgene[${c}]}.pydair fugu3.c${cgene[${c}]}.pydair \
                 -n Fugu1 Fugu2 Fugu3 \
                 -o fugustats_c${cgene[${c}]} \
                 --estimate-vdj-combination
done
```

Check the number of sequences and CDR3 that contains stop codons.

```
cgene=("m" "t")
for (( i = 1; i < 4; ++i ))
do
    for (( c = 0; c < ${#cgene[@]}; ++c ))
    do
        python ../bin/check_orf.py -i ./fugu${i}.c${cgene[${c}]}.pydair
    done
done
```

Create CDR3 sequence FASTA files.

```bash
cgene=("m" "t")         # cm, ct
for (( i = 1; i < 4; ++i ))
do
    for (( c = 0; c < ${#cgene[@]}; ++c ))
    do
        python ../bin/create_cdr3_fasta.py -i fugu${i}.c${cgene[${c}]}.pydair -f fugu${i}
    done
done
```

Perform rarefaction study for VDJ combinations for Cm group.

```
pwd
## ~/Desktop/fugu/work

Rscript --vanilla --slave ../bin/plot_rarefaction_vdjcombi.R \
          'fugustats_cm.Fugu1.vdj.freq.tsv,fugustats_cm.Fugu2.vdj.freq.tsv,fugustats_cm.Fugu3.vdj.freq.tsv' \
          'Fugu 1,Fugu 2,Fugu 3' Fig_Cm_VDJ_rarefaction.pdf pdf 10 6 100
```




## 6 Figures


Plot the frequencies of V, D, and J gene usages,
and the distribution of 3'V deletion, 5'J deletion, and 3'V5'J insertions.


```
pwd
## ~/Desktop/fugu/work

Rscript --vanilla --slave ../bin/plot_freq.R \
          'fugustats_cm.Fugu1,fugustats_cm.Fugu2,fugustats_cm.Fugu3' \
          'Fugu 1,Fugu 2,Fugu 3' TRUE Fig_Cm_VDJ_freq_hist.pdf pdf 14 4 300
Rscript --vanilla --slave ../bin/plot_freq.R \
          'fugustats_ct.Fugu1,fugustats_ct.Fugu2,fugustats_ct.Fugu3' \
          'Fugu 1,Fugu 2,Fugu 3' TRUE Fig_Ct_VDJ_freq_hist.pdf pdf 10 4 100

Rscript --vanilla --slave ../bin/plot_indelslen.R \
          fugustats_ct.all.v_del_len.freq.tsv TRUE Fig_Ct_Vdel.pdf pdf 5 4 10
Rscript --vanilla --slave ../bin/plot_indelslen.R \
          fugustats_cm.all.v_del_len.freq.tsv TRUE Fig_Cm_Vdel.pdf pdf 5 4 10

Rscript --vanilla --slave ../bin/plot_indelslen.R \
          fugustats_ct.all.j_del_len.freq.tsv TRUE Fig_Ct_Jdel.pdf pdf 8 4 20
Rscript --vanilla --slave ../bin/plot_indelslen.R \
          fugustats_cm.all.j_del_len.freq.tsv TRUE Fig_Cm_Jdel.pdf pdf 8 4 20

Rscript --vanilla --slave ../bin/plot_indelslen.R \
          fugustats_ct.all.vj_ins_len.freq.tsv TRUE Fig_Ct_VJins.pdf pdf 14 4 40
Rscript --vanilla --slave ../bin/plot_indelslen.R \
          fugustats_cm.all.vj_ins_len.freq.tsv TRUE Fig_Cm_VJins.pdf pdf 14 4 40

Rscript --vanilla --slave ../bin/plot_indelslen.R \
          fugustats_ct.all.cdr3_prot_len.freq.tsv TRUE Fig_Ct_CDR3aa.pdf pdf 8 4 20
Rscript --vanilla --slave ../bin/plot_indelslen.R \
          fugustats_cm.all.cdr3_prot_len.freq.tsv TRUE Fig_Cm_CDR3aa.pdf pdf 8 4 20

```


Plot 3D figures for visualizing VDJ combinations.
One should start R and change the variables to plot the 3D scatter plots for thee three fugu.

```bash
$R

>source('../bin/plot_vdj_3d.R')
```






## 7 CDR3 diversity study

CDR3 amino acid sequences of the three fugu are merged into one file for estimating CDR3 repertoire population sizes.
Then, CDR3 amino acid sequences were clustered by [CD-HIT](http://weizhongli-lab.org/cd-hit/). Output files from CD-HIT were adjusted by Python script, and estimation was performed by R script.


```bash
pwd
## ~/Desktop/fuug/work

cgene=("m" "t")


## merge fugu1, fugu2, fugu3 data into one file.
for ((c = 0; c < ${#cgene[@]}; ++c))
do
    cat fugu1.c${cgene[${c}]}.pydair.cdr3prot.fa > fuguall.c${cgene[${c}]}.cdr3prot.fa
    cat fugu2.c${cgene[${c}]}.pydair.cdr3prot.fa >> fuguall.c${cgene[${c}]}.cdr3prot.fa
    cat fugu3.c${cgene[${c}]}.pydair.cdr3prot.fa >> fuguall.c${cgene[${c}]}.cdr3prot.fa
done


## cluster CDR3 sequences using CD-HIT
for ((c = 0; c < ${#cgene[@]}; ++c))
do
    ../bin/cd-hit-v4.6.6-2016-0711/cd-hit \
             -i fuguall.c${cgene[${c}]}.cdr3prot.fa \
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


python ../bin/make_cdr3_clusters_tsv.py -f fuguall.cm.cdr3aa.fa          \
                                        -c cdhit_fuguall.cm.cls.fa.clstr \
                                        -o fuguall.cm.cdr3.venn.tsv
python ../bin/make_cdr3_clusters_tsv.py -f fuguall.ct.cdr3aa.fa          \
                                        -c cdhit_fuguall.ct.cls.fa.clstr \
                                        -o fuguall.ct.cdr3.venn.tsv                                        
```







