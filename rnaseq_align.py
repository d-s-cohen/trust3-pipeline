#!/usr/bin/env python

'''
RNA-Seq Alignment Workflow SLURM Script Generator
Generates and submits a SLURM script that takes a BAM file, converts it to FASTQ, uses STAR two-pass method via ICGC code, and produces an aligned BAM file
David Cohen, February 2019
https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
'''

import argparse
import os
import sys


parser = argparse.ArgumentParser(description="RNA-Seq Alignment Workflow SLURM Submission File Generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage='%(prog)s [options]', add_help=True)
required = parser.add_argument_group("required input parameters")
required.add_argument("--bamIn", default=None, help="Input unaligned BAM file", required=True)
optional = parser.add_argument_group("optional input parameters")
optional.add_argument("--out", help="String which output directory and file names are based upon. Default string is based on input file name.", required=False)
optional.add_argument("--workDir", default=os.path.abspath('.'), help="Work directory")

args = parser.parse_args()

# Define output name

if args.out is None:
	outName = os.path.basename(args.bamIn)
        outName = os.path.splitext(outName)[0]
else:
        outName = args.out

workDir = args.workDir

# workDir configuation

if workDir.endswith("/"):
        workDir = workDir[:-1]

job_file = workDir + "/Slurm/rnaseq_align_" + outName + ".sh"

# Write SLURM script

with open(job_file,"w") as fh:

        fh.writelines("#!/bin/bash\n")
        fh.writelines("\n")
        fh.writelines("#SBATCH --job-name=rnaseq_align\n")
        fh.writelines("#SBATCH --time=06:00:00\n")
        fh.writelines("#SBATCH --mem=40G\n")
        fh.writelines("#SBATCH --cpus-per-task=11\n")
        fh.writelines("#SBATCH --error=" + workDir + "/Slurm/" + outName + "_error.err\n")
        fh.writelines("#SBATCH --output=" + workDir + "/Slurm/" + outName + "_outfile.out\n")
        fh.writelines("\n")
        fh.writelines("###############################################################\n")
        fh.writelines("# RNA-Seq Alignment Workflow \n")
        fh.writelines("# - Takes BAM file, converts to FASTQ, uses STAR two-pass method via ICGC code, produces an aligned BAM\n")
        fh.writelines("# - David Cohen, February 2019\n")
        fh.writelines("# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/\n")
        fh.writelines("###############################################################\n")
        fh.writelines("\n")
        fh.writelines('[ -d "' + workDir + '/' + outName + '" ] && (echo "Directory ' + outName + ' already exists. Choose a Different output name or remove current directory."; scancel $SLURM_JOB_ID)\n')
        fh.writelines("\n")
        fh.writelines("mkdir " + workDir + "/" + outName + "\n")
        fh.writelines("\n")
        fh.writelines("source activate RNA-Seq_Alignment\n")
        fh.writelines("\n")
        fh.writelines("# BAM to FASTQ conversion via biobambam\n")
        fh.writelines("\n")
        fh.writelines("bamtofastq \\\n")
        fh.writelines("collate=1 \\\n")
        fh.writelines("exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \\\n")
        fh.writelines("filename=" + args.bamIn + " \\\n")
        fh.writelines("outputdir=" + workDir + "/" + outName + " \\\n")
        fh.writelines("outputperreadgroup=1 \\\n")
        fh.writelines("gz=1 \\\n")
        fh.writelines("outputperreadgroupsuffixF=_1.fastq.gz \\\n")
        fh.writelines("outputperreadgroupsuffixF2=_2.fastq.gz \\\n")
        fh.writelines("tryoq=1 \\\n")
        fh.writelines("T=" + workDir + "/" + outName + "/" + "collation.TEMP \\\n")
        fh.writelines("inputformat=bam\n")
        fh.writelines("\n")
        fh.writelines("# Prepare FASTQ for star_align.py\n")
        fh.writelines("\n")
        fh.writelines("tar cvf " + workDir + "/" + outName + "/" + outName + "_fastq.tar \\\n")
        fh.writelines(workDir + "/" + outName + "/" + outName + "_1.fastq.gz \\\n")
        fh.writelines(workDir + "/" + outName + "/" + outName + "_2.fastq.gz \\\n")
        fh.writelines("--remove-files\n")
        fh.writelines("\n")
        fh.writelines("# Run Star Align \n")
        fh.writelines("\n")
        fh.writelines("python " + workDir + "/Resources/icgc_rnaseq_align/star_align.py \\\n")
        fh.writelines("--genomeDir " + workDir + "/Resources/Index/star_genome_d1_vd1_gtfv22 \\\n")
        fh.writelines("--tarFileIn " + workDir + "/" + outName + "/" + outName + "_fastq.tar \\\n")
        fh.writelines("--workDir " + workDir + "/" + outName + " \\\n")
        fh.writelines("--genomeFastaFiles " + workDir + "/Resources/Genome/GRCh38.d1.vd1.fa \\\n")
        fh.writelines("--out " + workDir + "/" + outName + "/" + outName + "_aligned.bam \\\n")
        fh.writelines("--runThreadN 8 \\\n")
        fh.writelines("--outFilterMultimapScoreRange 1 \\\n")
        fh.writelines("--outFilterMultimapNmax 20 \\\n")
        fh.writelines("--outFilterMismatchNmax 10 \\\n")
        fh.writelines("--alignIntronMax 500000 \\\n")
        fh.writelines("--alignMatesGapMax 1000000 \\\n")
        fh.writelines("--sjdbScore 2 \\\n")
        fh.writelines("--limitBAMsortRAM 0 \\\n")
        fh.writelines("--alignSJDBoverhangMin 1 \\\n")
        fh.writelines("--genomeLoad NoSharedMemory \\\n")
        fh.writelines("--outFilterMatchNminOverLread 0.33 \\\n")
        fh.writelines("--outFilterScoreMinOverLread 0.33 \\\n")
        fh.writelines("--twopass1readsN -1 \\\n")
        fh.writelines("--sjdbOverhang 100 \\\n")
        fh.writelines("--outSAMstrandField intronMotif \\\n")
        fh.writelines("--outSAMunmapped Within\n")
        fh.writelines("\n")
        fh.writelines("source deactivate\n")
        fh.writelines("\n"

fh.close()

# Submit SLURM script

os.system("sbatch " + job_file)


