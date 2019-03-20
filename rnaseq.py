#!/usr/bin/env python

'''
RNA-Seq Alignment and TRUST Analysis SLURM Script Generator
Generates and submits a SLURM script that takes a BAM file, converts it to FASTQ, uses STAR two-pass method via ICGC code, produces an aligned BAM and BAM index file, and runs analysis with TRUST
Also can take a directory of BAM files and perform the process on each BAM using a SLURM array
David Cohen, February - March 2019
https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
https://www.nature.com/articles/ng.3820
'''

import argparse
import os
import sys
import errno

parser = argparse.ArgumentParser(description="RNA-Seq Alignment and TRUST Analysis SLURM Script Generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage='%(prog)s [options]', add_help=True)
required = parser.add_argument_group("required input parameters")
required.add_argument("-i", "--bamIn", default=None, help="Either a directory containing unaligned BAM files or a single BAM file", required=True)
optional = parser.add_argument_group("optional input parameters")
optional.add_argument("-o", "--out", help="String which SLURM file names are based upon in case of input directory. String which temporary directory and all output file names are based upon in case of single input file. Default string is based on input name.", required=False)
optional.add_argument("-w", "--workDir", default=os.path.abspath('.'), help="Work directory")
optional.add_argument("-e", "--endsWith", default="", help="For directory input, scan for files ending in .bam, preceded by this string, in order to specifiy a subset of .bam files.")
optional.add_argument('-s', action='store_true', help="If selected, script will be generated but not submitted to slurm")

args = parser.parse_args()

if  os.path.isfile(args.bamIn):

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

	job_file = workDir + "/slurm/" + outName + ".sh"

	# Write SLURM script

	with open(job_file,"w") as fh:

		fh.writelines("#!/bin/bash\n")
		fh.writelines("\n")
		fh.writelines("#SBATCH --job-name=rnaseq\n")
		fh.writelines("#SBATCH --time=06:00:00\n")
		fh.writelines("#SBATCH --mem=45G\n")
		fh.writelines("#SBATCH --cpus-per-task=11\n")
		fh.writelines("#SBATCH --error=" + workDir + "/slurm/" + outName + ".err\n")
		fh.writelines("#SBATCH --output=" + workDir + "/slurm/" + outName + ".out\n")
		fh.writelines("\n")
		fh.writelines("###############################################################\n")
		fh.writelines("# RNA-Seq Alignment and TRUST Analysis \n")
		fh.writelines("# - Takes BAM file, converts to FASTQ, uses STAR two-pass method via ICGC code, produces an aligned BAM and BAM index file, and runs analysis with TRUST\n")
		fh.writelines("# - David Cohen, February - March 2019\n")
		fh.writelines("# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/\n")
		fh.writelines("# https://www.nature.com/articles/ng.3820\n")
		fh.writelines("###############################################################\n")
		fh.writelines("\n")
		fh.writelines('[ -d "' + workDir + '/' + outName + '" ] && (echo "Directory ' + outName + ' already exists. Choose a Different output name or remove current directory."; scancel $SLURM_JOB_ID)\n')
		fh.writelines("\n")
		fh.writelines("mkdir " + workDir + "/" + outName + "\n")
		fh.writelines("mkdir " + workDir + "/" + outName + "/fastq\n")
		fh.writelines("\n")
		fh.writelines("source activate RNA-Seq_Alignment\n")
		fh.writelines("\n")
		fh.writelines("# BAM to FASTQ conversion via biobambam\n")
		fh.writelines("\n")
		fh.writelines("bamtofastq \\\n")
		fh.writelines("collate=1 \\\n")
		fh.writelines("exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \\\n")
		fh.writelines("filename=" + args.bamIn + " \\\n")
		fh.writelines("outputdir=" + workDir + "/" + outName + "/fastq \\\n")
		fh.writelines("outputperreadgroup=1 \\\n")
		fh.writelines("gz=1 \\\n")
		fh.writelines("outputperreadgroupsuffixF=_1.fastq.gz \\\n")
		fh.writelines("outputperreadgroupsuffixF2=_2.fastq.gz \\\n")
		fh.writelines("tryoq=1 \\\n")
		fh.writelines("T=" + workDir + "/" + outName + "/" + "collation.TEMP \\\n")
		fh.writelines("inputformat=bam\n")
		fh.writelines("\n")
		fh.writelines("# Run Star Align \n")
		fh.writelines("\n")
		fh.writelines("python " + workDir + "/resources/icgc_rnaseq_align/star_align.py \\\n")
		fh.writelines("--genomeDir " + workDir + "/resources/index/star_genome_d1_vd1_gtfv22 \\\n")
		fh.writelines("--fastqDir " + workDir + "/" + outName + "/fastq \\\n")
		fh.writelines("--workDir " + workDir + "/" + outName + " \\\n")
		fh.writelines("--genomeFastaFiles " + workDir + "/resources/genome/GRCh38.d1.vd1.fa \\\n")
		fh.writelines("--annotation " + workDir + "/resources/annotation/gencode.v22.annotation.gtf \\\n")
		fh.writelines("--out " + workDir + "/" + outName + "/" + outName + ".bam \\\n")
		fh.writelines("--runThreadN 8 \\\n")
		fh.writelines("--outFilterMultimapScoreRange 1 \\\n")
		fh.writelines("--outFilterMultimapNmax 20 \\\n")
		fh.writelines("--outFilterMismatchNmax 10 \\\n")
		fh.writelines("--alignIntronMax 500000 \\\n")
		fh.writelines("--alignMatesGapMax 1000000 \\\n")
		fh.writelines("--sjdbScore 2 \\\n")
		fh.writelines("--limitBAMsortRAM 44000000000 \\\n")
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
		fh.writelines("samtools index " +  workDir + "/" + outName + "/" + outName + ".bam " +  workDir + "/" + outName + "/" + outName + ".bam.bai \n")
		fh.writelines("\n")
		fh.writelines("# Sort everything \n")
		fh.writelines("\n")
		fh.writelines("mv  -t " + workDir + "/bam_aligned/ " + workDir + "/" + outName + "/" + outName + ".bam " +  workDir + "/" + outName + "/" + outName + ".bam.bai \n" )
		fh.writelines("mv  -t " + workDir + "/fastq/ " + workDir + "/" + outName + "/fastq/" + outName + "_1.fastq.gz " + workDir + "/" + outName + "/fastq/" + outName + "_2.fastq.gz \n")
		fh.writelines("mv " + workDir + "/" + outName + "/Log.final.out " + workDir + "/log/" + outName + ".out \n")
		fh.writelines("mv " + workDir + "/" + outName + "/Log_1st_pass.final.out " + workDir + "/log/" + outName + "_1st_pass.out \n")
		fh.writelines("rm -r " + workDir + "/" + outName + "\n")
		fh.writelines("\n")
		fh.writelines("# Run TRUST \n")
		fh.writelines("\n")
		fh.writelines("trust -f " + workDir + "/bam_aligned/" + outName + ".bam -g hg38 -E -o " + workDir + "/trust/ \n")
		fh.writelines("trust -f " + workDir + "/bam_aligned/" + outName + ".bam -g hg38 -B -o " + workDir + "/trust/ \n")
		fh.writelines("trust -f " + workDir + "/bam_aligned/" + outName + ".bam -g hg38 -B -L -o " + workDir + "/trust/ \n")
		fh.writelines("\n")
		fh.writelines("source deactivate\n")
		fh.writelines("\n")

	fh.close()

	# Submit SLURM script

        if not args.s:
		os.system("sbatch " + job_file)

elif os.path.isdir(args.bamIn):

	bamNameList = []

	for file in os.listdir(args.bamIn):
		if file.endswith(args.endsWith + ".bam"):
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

	# Directory structure setup

	try:
		os.makedirs(workDir + "/" +  outName)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	try:
		os.makedirs(workDir + "/" +  outName + "/bam")
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	try:
		os.makedirs(workDir + "/" +  outName + "/fastq")
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	try:
		os.makedirs(workDir + "/" +  outName + "/trust")
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	try:
		os.makedirs(workDir + "/" +  outName + "/slurm")
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	try:
		os.makedirs(workDir + "/" +  outName + "/log")
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	try:
		os.makedirs(workDir + "/" +  outName + "/TEMP")
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	try:
		os.makedirs(workDir + "/" +  outName + "/htseq")
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	job_file = workDir + "/" + outName + "/slurm/" + outName + ".sh"

	# Write SLURM script

	with open(job_file,"w") as fh:

		fh.writelines("#!/bin/bash\n")
		fh.writelines("\n")
		fh.writelines("#SBATCH --job-name=rnaseq\n")
		fh.writelines("#SBATCH --time=10:00:00\n")
		fh.writelines("#SBATCH --mem=45G\n")
		fh.writelines("#SBATCH --cpus-per-task=11\n")
		fh.writelines("#SBATCH --error=" + workDir + "/" + outName + "/" + "slurm/%A_%a.err\n")
		fh.writelines("#SBATCH --output=" + workDir + "/" + outName + "/" + "slurm/%A_%a.out\n")
		fh.writelines("\n")
		fh.writelines("###############################################################\n")
		fh.writelines("# RNA-Seq Alignment and TRUST Analysis \n")
		fh.writelines("# Takes a directory of BAM files and performs the following process on each BAM using a SLURM array: \n")
		fh.writelines("# - Takes BAM file, converts to FASTQ, uses STAR two-pass method via ICGC code, produces an aligned BAM and BAM index file, and runs analysis with TRUST\n")
		fh.writelines("# - David Cohen, February - March 2019\n")
		fh.writelines("# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/\n")
		fh.writelines("# https://www.nature.com/articles/ng.3820\n")
		fh.writelines("###############################################################\n")
		fh.writelines("\n")
		fh.writelines('bamNameList=("' + '" "'.join(bamNameList) + '") \n')
                fh.writelines("# bash array length is " + str(arrayLength) + " \n")
		fh.writelines("\n")
		fh.writelines('[ -d "' + workDir + '/' + outName + '/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}" ] && (echo "Directory ${bamNameList[$SLURM_ARRAY_TASK_ID]} already exists. Choose a Different output name or remove current directory."; scancel $SLURM_JOB_ID)\n')
		fh.writelines("\n")
		fh.writelines("mkdir " + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}\n")
		fh.writelines("mkdir " + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq\n")
		fh.writelines("\n")
		fh.writelines("source activate RNA-Seq_Alignment\n")
		fh.writelines("\n")
		fh.writelines("# BAM to FASTQ conversion via biobambam\n")
		fh.writelines("\n")
		fh.writelines("bamtofastq \\\n")
                fh.writelines("F=" + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq/${bamNameList[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz \\\n")
                fh.writelines("F2=" + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq/${bamNameList[$SLURM_ARRAY_TASK_ID]}_2.fastq.gz \\\n")
		fh.writelines("collate=1 \\\n")
		fh.writelines("exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \\\n")
		fh.writelines("filename=" + bamDir + "/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam \\\n")
		fh.writelines("outputperreadgroup=0 \\\n")
		fh.writelines("gz=1 \\\n")
		fh.writelines("tryoq=1 \\\n")
		fh.writelines("T=" + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/" + "collation.TEMP \\\n")
		fh.writelines("inputformat=bam\n")
		fh.writelines("\n")
		fh.writelines("# Run Star Align \n")
		fh.writelines("\n")
		fh.writelines("python " + workDir + "/resources/icgc_rnaseq_align/star_align.py \\\n")
		fh.writelines("--genomeDir " + workDir + "/resources/index/star_genome_d1_vd1_gtfv22 \\\n")
		fh.writelines("--fastqDir " + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq \\\n")
		fh.writelines("--workDir " + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]} \\\n")
		fh.writelines("--genomeFastaFiles " + workDir + "/resources/genome/GRCh38.d1.vd1.fa \\\n")
		fh.writelines("--annotation " + workDir + "/resources/annotation/gencode.v22.annotation.gtf \\\n")
		fh.writelines("--out " + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam \\\n")
		fh.writelines("--runThreadN 8 \\\n")
		fh.writelines("--outFilterMultimapScoreRange 1 \\\n")
		fh.writelines("--outFilterMultimapNmax 20 \\\n")
		fh.writelines("--outFilterMismatchNmax 10 \\\n")
		fh.writelines("--alignIntronMax 500000 \\\n")
		fh.writelines("--alignMatesGapMax 1000000 \\\n")
		fh.writelines("--sjdbScore 2 \\\n")
		fh.writelines("--limitBAMsortRAM 44000000000 \\\n")
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
		fh.writelines("samtools index " +  workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam \\\n")
		fh.writelines(workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam.bai \n")
		fh.writelines("\n")
		fh.writelines("# Sort everything \n")
		fh.writelines("\n")
		fh.writelines("mv  -t " + workDir + "/" + outName + "/bam/ \\\n")
		fh.writelines(workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam \\\n")
		fh.writelines(workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam.bai \n")
		fh.writelines("\n")
		fh.writelines("mv  -t " + workDir + "/" + outName + "/fastq/ \\\n")
		fh.writelines(workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/fastq/* \\\n")
		fh.writelines("\n")
		fh.writelines("mv " + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/Log.final.out \\\n")
		fh.writelines(workDir + "/" + outName + "/log/${bamNameList[$SLURM_ARRAY_TASK_ID]}.out \n")
		fh.writelines("\n")
		fh.writelines("mv " + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}/Log_1st_pass.final.out \\\n")
		fh.writelines(workDir + "/" + outName + "/log/${bamNameList[$SLURM_ARRAY_TASK_ID]}_1st_pass.out \n")
		fh.writelines("\n")
		fh.writelines("rm -r " + workDir + "/" + outName + "/TEMP/${bamNameList[$SLURM_ARRAY_TASK_ID]}\n")
		fh.writelines("\n")
		fh.writelines("# Run TRUST \n")
		fh.writelines("\n")
		fh.writelines("trust -f " + workDir + "/" + outName + "/bam/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam -g hg38 -E -o " + workDir + "/" + outName + "/trust/ \n")
		fh.writelines("trust -f " + workDir + "/" + outName + "/bam/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam -g hg38 -B -o " + workDir + "/" + outName + "/trust/ \n")
		fh.writelines("trust -f " + workDir + "/" + outName + "/bam/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam -g hg38 -B -L -o " + workDir + "/" + outName + "/trust/ \n")
		fh.writelines("\n")
#		fh.writelines("htseq-count -m intersection-nonempty -i gene_id -r pos -s no -f bam " + workDir + "/" + outName + "/bam/${bamNameList[$SLURM_ARRAY_TASK_ID]}.bam " + workDir + "/resources/annotation/gencode.v22.annotation.gtf > " + workDir + "/" + outName + "/htseq/${bamNameList[$SLURM_ARRAY_TASK_ID]}.txt \n")		
		fh.writelines("\n")
		fh.writelines("source deactivate\n")
		fh.writelines("\n")

	fh.close()

	# Submit SLURM script

	if not args.s:
		os.system("sbatch --array=0-" + str(arrayLength) + " " + job_file)

else:

	print("Input is neither a directory or a file")

