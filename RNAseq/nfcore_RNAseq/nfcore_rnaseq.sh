#!/bin/bash
#SBATCH --job-name=rna
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --time=3-00:00:0
#SBATCH --mem=30G
#SBATCH --partition=ncpu
#SBATCH --mail-type="END,FAIL"
#SBATCH --mail-user=nicolas.veland@crick.ac.uk
#SBATCH --output=run_rnaseq.o
#SBATCH --error=run_rnaseq.e

# many of the changes you need to make to this file are the same as those 
# changes made to the '1.fetchngs-job-template-rnaseq.sh' file. Please refer to
# the comments in '1.fetchngs-job-template-rnaseq.sh' for tips on modifying
# this file.

# Load software modules - check for newer versions using 'ml spider Nextflow'
ml purge
ml Nextflow/25.10.0
ml Singularity/3.11.3

# Create a 'nextflow' folder in your working area. Run 'cd nextflow' and 'pwd -P'
# Copy the result to replace the directory path below
export NXF_HOME=/nemo/lab/haydaya/home/users/velandn/working/TCR_RNAseq_1/rnaseq/nextflow/


# NXF_WORK will be created in the directory you execute this script from
export NXF_WORK=/nemo/lab/haydaya/home/users/velandn/working/TCR_RNAseq_1/rnaseq/work/


# You can leave this line as is - we are using the BABS cache of singularity images
export NXF_SINGULARITY_LIBRARYDIR=/flask/apps/containers/all-singularity-images/

# You should create a folder in your working area to store any additional images
# that are not found in the BABS library folder above
export NXF_SINGULARITY_CACHEDIR=/nemo/lab/haydaya/home/users/velandn/working/TCR_RNAseq_1/rnaseq/cache/


# run the pipeline
nextflow run main.nf \
    --input samplesheet.csv  \
    --outdir rnaseq-results_1/ \
    --aligner star_rsem \
    --fasta Mus_musculus.GRCm39.dna_sm.primary_assembly.fa \
    --gtf Mus_musculus.GRCm39.115.gtf \
    -c nextflow.config \
    --skip_biotype_qc \
    --save_reference \
    -profile crick \
    -resume


# Once to have made your changes to this file, save it and submit the job with:
#   sbatch 2.rnaseq-job-template.sh 
