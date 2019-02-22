## RNA-Seq Alignment Workflow SLURM Script Generator

#### David Cohen - February 2019

Generates and submits a SLURM script that takes a BAM file, converts it to FASTQ, uses STAR two-pass method via ICGC code, and produces an aligned BAM file

Conducts work in an individual generated temporary directory for each job and then distributes output BAM, FASTQ, and other files to appropriate directories  

##### Usage:

rnaseq_align.py [options]

```
required input parameters:
  --bamIn BAMIN      Input unaligned BAM file (default: None)

optional input parameters:
  --out OUT          String which temporary directory and output file names are based
                     upon. Default string is based on input file name
                     (default: None)
  --workDir WORKDIR  Work directory (default: ./)
```

##### Requirements:

Attached python environment RNA-Seq_Alignment:

```
conda env create -f resources/environment.yml
```

File structure in working directory:

* resources/icgc_rnaseq_align/star_align.py modified from <https://github.com/akahles/icgc_rnaseq_align>

* resources/genome/GRCh38.d1.vd1.fa from <https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files> (GRCh38.d1.vd1 Reference Sequence)

* resources/index/star_genome_d1_vd1_gtfv22/ directory and contents from <https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files> (GDC.h38.d1.vd1 STAR2 Index Files)

* resources/annotation/gencode.v22.annotation.gtf from <https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files> (GDC.h38 GENCODE v22 GTF)

##### Reference: 

<https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>

