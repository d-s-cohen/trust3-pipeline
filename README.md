## RNA-Seq Alignment Workflow Slurm Script Generator and Submitter

#### David Cohen - February 2019

Generates and submits a Slurm script that takes a BAM file, converts it to FASTQ, uses STAR two-pass method via ICGC code, and produces an aligned BAM file

##### Usage:

rnaseq_align.py [options]

```
required input parameters:
  --bamIn BAMIN      Input unaligned BAM file (default: None)

optional input parameters:
  --out OUT          String which output directory and file names are based
                     upon. Default string is based on input file name
                     (default: None)
  --workDir WORKDIR  Work directory (default: ./)
```

##### Requirements:

Included conda environment RNA-Seq_Alignment:

```
conda env create -f Resources/environment.yml
```

File structure in working directory:

* Resources/icgc_rnaseq_align/star_align.py which is modified from <https://github.com/akahles/icgc_rnaseq_align> and already included in this repository

* Resources/Genome/GRCh38.d1.vd1.fa from <https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files> (GRCh38.d1.vd1 Reference Sequence)

* Resources/Index/star_genome_d1_vd1_gtfv22/ directory and contents from <https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files> (GDC.h38.d1.vd1 STAR2 Index Files)

SAMtools

##### Reference: 

<https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>

