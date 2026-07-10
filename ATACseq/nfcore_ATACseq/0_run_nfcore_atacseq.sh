
# nf-core/atacseq v2.0 run (mouse GRCm38, MACS2 broad peaks)


BASE=/path/to/project
REF=${BASE}/references
EMAIL=email@example.com

export TMPDIR=${BASE}/tmp
export NXF_TEMP=${BASE}/tmp
mkdir -p ${BASE}/singularity_cache ${TMPDIR}

cd ${BASE}

nextflow run nf-core/atacseq -r 2.0 \
  -profile singularity \
  -c nextflow.config \
  --input ${BASE}/download/samplesheet_nfcore_skin_clean.csv \
  --fasta ${REF}/Mus_musculus.GRCm38.dna_sm.toplevel.fa \
  --gtf ${REF}/Mus_musculus.GRCm38.95.gtf \
  --bwa_index ${REF}/ \
  --gene_bed ${REF}/Mus_musculus.GRCm38.95.bed \
  --blacklist ${REF}/mm10-blacklist.GRCm38.bed \
  --macs_gsize 2600000000 \
  --mito_name MT \
  --deseq2_vst \
  --min_reps_consensus 3 \
  --outdir results_atacseq \
  -work-dir werk/ \
  --email_on_fail ${EMAIL} \
  -resume
