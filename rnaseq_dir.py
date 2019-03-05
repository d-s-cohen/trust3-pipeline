#!/usr/bin/env python

'''
RNA-Seq Alignment and Analysis Workflow SLURM Script Generator
Takes a directory of BAM files and performs the following process on each BAM using a SLURM array:
Generates and submits a SLURM script that takes a BAM file, converts it to FASTQ, uses STAR two-pass method via ICGC code, produces an aligned BAM and BAM index file, and runs analysis with TRUST
David Cohen, February 2019
https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
'''

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="RNA-Seq Alignment Workflow SLURM Submission File Generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage='%(prog)s [options]', add_help=True)
required = parser.add_argument_group("required input parameters")
required.add_argument("-i", "--bamIn", default=None, help="Input directory containing unaligned BAM files", required=True)
optional = parser.add_argument_group("optional input parameters")
optional.add_argument("-o", "--out", help="String which SLURM file names are based upon.. Default string is based on input directory name.", required=False)
optional.add_argument("-w", "--workDir", default=os.path.abspath('.'), help="Work directory")

args = parser.parse_args()

# Define output name

bamNameList = []

for file in os.listdir(args.bamIn):
    if file.endswith(".bam"):
        bamNameList.append(file[:-4])

arrayLength = len(bamNameList) - 1

if args.out is None:
        outName = os.path.basename(os.path.normpath(args.bamIn))
else:
        outName = args.out

# workDir configuation

workDir = args.workDir

if workDir.endswith("/"):
        workDir = workDir[:-1]

# bamDir configuration

bamDir = os.path.abspath(args.bamIn)

if bamDir.endswith("/"):
        bamDir = bamDir[:-1]

job_file = workDir + "/slurm/" + outName + ".sh"

# Write SLURM script

with open(job_file,"w") as fh:

        fh.writelines("#!/bin/bash\n")
        fh.writelines("\n")
        fh.writelines("#SBATCH --job-name=rnaseq_align\n")
        fh.writelines("#SBATCH --time=06:00:00\n")
        fh.writelines("#SBATCH --mem=40G\n")
        fh.writelines("#SBATCH --cpus-per-task=11\n")
        fh.writelines("#SBATCH --error=" + workDir + "/slurm/" + outName + "_%A_%a.err\n")
        fh.writelines("#SBATCH --output=" + workDir + "/slurm/" + outName + "_%A_%a.out\n")
        fh.writelines("\n")
        fh.writelines("###############################################################\n")
        fh.writelines("# RNA-Seq Alignment and Analysis Workflow \n")
        fh.writelines("# Takes a directory of BAM files and performs the following process on each BAM using a SLURM array: \n")
        fh.writelines("# - Takes BAM file, converts to FASTQ, uses STAR two-pass method via ICGC code, produces an aligned BAM and BAM index file, and runs analysis with TRUST\n")
        fh.writelines("# - David Cohen, February 2019\n")
        fh.writelines("# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/\n")
        fh.writelines("###############################################################\n")
        fh.writelines("\n")
        fh.writelines('bamNameList=("' + '" "'.join(bamNameList) + '") \n')
        fh.writelines("\n")
        fh.writelines('[ -d "' + workDir + '/${bamNameList[$SLURM_ARRAY_TASK_ID]}" ] && (echo "Directory ${bamNameList[$SLURM_ARRAY_TASK_ID]} already exists. Choose a Different output name or remove current directory."; scancel $SLURM_JOB_ID)\n')
        fh.writelines("\n")
        fh.writelines("mkdir " + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}\n")
        fh.writelines("mkdir " + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq\n")
        fh.writelines("\n")
        fh.writelines("source activate RNA-Seq_Alignment\n")
        fh.writelines("\n")
        fh.writelines("# BAM to FASTQ conversion via biobambam\n")
        fh.writelines("\n")
        fh.writelines("bamtofastq \\\n")
        fh.writelines("collate=1 \\\n")
        fh.writelines("exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \\\n")
        fh.writelines("filename=" + bamDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam \\\n")
        fh.writelines("outputdir=" + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq \\\n")
        fh.writelines("outputperreadgroup=1 \\\n")
        fh.writelines("gz=1 \\\n")
        fh.writelines("outputperreadgroupsuffixF=_1.fastq.gz \\\n")
        fh.writelines("outputperreadgroupsuffixF2=_2.fastq.gz \\\n")
        fh.writelines("tryoq=1 \\\n")
        fh.writelines("T=" + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/" + "collation.TEMP \\\n")
        fh.writelines("inputformat=bam\n")
        fh.writelines("\n")
        fh.writelines("# Run Star Align \n")
        fh.writelines("\n")
        fh.writelines("python " + workDir + "/resources/icgc_rnaseq_align/star_align.py \\\n")
        fh.writelines("--genomeDir " + workDir + "/resources/index/star_genome_d1_vd1_gtfv22 \\\n")
        fh.writelines("--fastqDir " + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq \\\n")
        fh.writelines("--workDir " + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]} \\\n")
        fh.writelines("--genomeFastaFiles " + workDir + "/resources/genome/GRCh38.d1.vd1.fa \\\n")
        fh.writelines("--annotation " + workDir + "/resources/annotation/gencode.v22.annotation.gtf \\\n")
        fh.writelines("--out " + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam \\\n")
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
        fh.writelines("# Generate BAM index file\n")
        fh.writelines("\n")
        fh.writelines("samtools index " +  workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam " +  workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam.bai \n")
        fh.writelines("\n")
        fh.writelines("# Sort everything \n")
        fh.writelines("\n")
        fh.writelines("mv  -t " + workDir + "/bam_aligned/ \\\n")
        fh.writelines(workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam \\\n")
        fh.writelines(workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam.bai \n")
        fh.writelines("\n")
        fh.writelines("mv  -t " + workDir + "/fastq/ \\\n")
        fh.writelines(workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq/${bamNameList[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz \\\n")
        fh.writelines(workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq/${bamNameList[$SLURM_ARRAY_TASK_ID]}_2.fastq.gz \n")
        fh.writelines("\n")
        fh.writelines("mv " + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/Log.final.out \\\n")
        fh.writelines(workDir + "/log/${bamNameList[$SLURM_ARRAY_TASK_ID]}.out \n")
        fh.writelines("\n")
        fh.writelines("mv " + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}/Log_1st_pass.final.out \\\n")
        fh.writelines(workDir + "/log/${bamNameList[$SLURM_ARRAY_TASK_ID]}_1st_pass.out \n")
        fh.writelines("\n")
        fh.writelines("rm -r " + workDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}\n")
        fh.writelines("\n")
        fh.writelines("# Run TRUST \n")
        fh.writelines("\n")
        fh.writelines("trust -f " + workDir + "/bam_aligned/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam -g hg38 -o " + workDir + "/trust/ \n")
        fh.writelines("\n")
        fh.writelines("source deactivate\n")
        fh.writelines("\n")

fh.close()

# Submit SLURM script

os.system("sbatch --array=0-" + str(arrayLength) + " " + job_file)
