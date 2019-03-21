## RNA-Seq Alignment and TRUST Analysis Pipeline

#### David Cohen | February - March 2019

Generates and submits a SLURM script that takes an unaligned BAM file, converts it to FASTQ, uses STAR two-pass method via ICGC code, produces an aligned BAM file, then runs TRUST for analysis.

Conducts conversion and alignment in an individual generated temporary directory for each job and then distributes output BAM, FASTQ, and other files to appropriate directories. Then runs TRUST.

Capable of processing individual BAM files as input or directories of BAMs.  

##### Usage:

rnaseq.py [options]

```
required input parameters:
  -i BAMIN, --bamIn BAMIN
                        Either a directory containing unaligned BAM files or a
                        single BAM file (default: None)

optional input parameters:
  -o OUT, --out OUT     String which SLURM file names are based upon in case
                        of input directory. String which temporary directory
                        and all output file names are based upon in case of
                        single input file. Default string is based on input
                        name. (default: None)
  -w WORKDIR, --workDir WORKDIR
                        Work directory (default: ./)
  -e EXT, --ext EXT     For directory input, scan for files ending in this
                        string (default: .bam)
  -p EXTPRE, --extPre EXTPRE
                        For directory input, this is an optional string to
                        precede the extension, in order to specifiy a subset
                        of those files. (default: )
  -s                    If selected, script will be generated but not
                        submitted to slurm. Useful for modifying the script
                        before submission. (default: False)
```

##### Requirements:

Create python environment RNA-Seq_Alignment from attached environment file:

```
conda env create -f resources/environment.yml
```

Then install TRUST and htseq in environment:

Download the latest version of TRUST from <https://bitbucket.org/liulab/trust>

```
tar xvzf trust-*.tar.gz
cd trust
source activate RNA-Seq_Alignment
python setup.py install
conda install htseq
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
