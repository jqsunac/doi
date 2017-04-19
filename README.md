# Timecourse RNA-Seq data analysis of C. insueta


## QC

```
cd qc/
bash qc.sh
```

## Assembly

To assembly A-genome and R-genome:

```
cd assembly/
bash refassemble_rnaseq.sh camara
bash refassemble_rnaseq.sh crivularis
```

To calculate statistics:


```
bash calc_stats.sh camara
bash calc_stats.sh crivularis

```


## Mapping


To mapping, read classification, counting:

```
cd rnaseq_processing/
bash proc_main.sh

```




