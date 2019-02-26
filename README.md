## RNA-Seq Alignment and TRUST Analysis Workflow

#### David Cohen - February 2019

Generates and submits a SLURM script that takes an unaligned BAM file, converts it to FASTQ, uses STAR two-pass method via ICGC code, produces an aligned BAM file, then runs TRUST for analysis.

Conducts work in an individual generated temporary directory for each job and then distributes output BAM, FASTQ, and other files to appropriate directories  

##### Usage:

rnaseq_align.py [options]

```
required input parameters:
  -i BAMIN, --bamIn BAMIN
                        Input unaligned BAM file (default: None)

optional input parameters:
  -o OUT, --out OUT     String which temporary directory and output file names
                        are based upon. Default string is based on input file
                        name. (default: None)
  -w WORKDIR, --workDir WORKDIR
                        Work directory (default: ./)
```

##### Requirements:

Create python environment RNA-Seq_Alignment from attached environment file:

```
conda env create -f resources/environment.yml
```

Then install TRUST in environment:

Download the latest version of TRUST from <https://bitbucket.org/liulab/trust>

```
tar xvzf trust-*.tar.gz
cd trust
source activate RNA-Seq_Alignment
python setup.py install
source deactivate
```

File structure in working directory:

* resources/icgc_rnaseq_align/star_align.py modified from <https://github.com/akahles/icgc_rnaseq_align>

* resources/genome/GRCh38.d1.vd1.fa from <https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files> (GRCh38.d1.vd1 Reference Sequence)

* resources/index/star_genome_d1_vd1_gtfv22/ directory and contents from <https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files> (GDC.h38.d1.vd1 STAR2 Index Files)

* resources/annotation/gencode.v22.annotation.gtf from <https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files> (GDC.h38 GENCODE v22 GTF)

SAMtools - <http://www.htslib.org/download/>

##### Reference: 

<https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>

<https://www.nature.com/articles/ng.3820>
