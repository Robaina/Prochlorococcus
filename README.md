This repository contains the python code employed to generate the environmental normalized count data shown in the accompanying publication.


## Installing conda environment

1. Download [diffexpr](https://github.com/wckdouglas/diffexpr) python repository (runs [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in Python) from GitHub

2. Install and activate conda environment:

```bash
conda create env -f metat2genome.yml

conda activate metat2genome
```

3. Install _diffexpr_:

```bash
cd diffexpr
Rscript setup.R
python setup.py install
```

## Order and function of scripts:

1. ```filter_fastq_by_names.py```: filter reads in metagenomic fastq files by blastx (against MIT9301) results
2. ```align_fastq.py```: aligns metagenomic fastq files to MIT9301 reference genome, then extracts alignments containing MD tags from SAM files (MD tags are required to compute percent identity later on)
3. ```get_filtered_counts.py```: filters alignments by read length, percent of matched sequence and percent identity.
4. ```get_normalized_counts.py```: deseq2 normalization of counts and sorting of samples based on environmental factors
