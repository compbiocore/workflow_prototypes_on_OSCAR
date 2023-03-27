### 1. Pull Down the Repo

If you haven't already pulled down the repo

```bash 
cd $HOME
git clone https://github.com/compbiocore/workflows_on_OSCAR.git
```

### 2. Modify the CSV File
Edit your CSV file

```commandline
## Fix the CSV file
sample_id,read1,read2
S1,/gpfs/data/cbc/sanghyun_lee/01.RawData/S1/S1_1.fq.gz,/gpfs/data/cbc/sanghyun_lee/01.RawData/S1/S1_2.fq.gz
S2,/gpfs/data/cbc/sanghyun_lee/01.RawData/S2/S2_1.fq.gz,/gpfs/data/cbc/sanghyun_lee/01.RawData/S2/S2_2.fq.gz
```

### 3. Run the RNASeq

```commandline
cd /gpfs/data/cbc/sanghyun_lee/
nextflow_start

nextflow run /gpfs/data/cbc/pcao5/workflow_prototypes_on_OSCAR/workflows/rnaseq.nf \
--samplesheet samplesheet_sanghyn_test.csv --out_dir out \
--reference_genome Mus_musculus.GRCm38.fa --gtf gencode.vM25.annotation.gtf
```

### 4. Run the HOMER Pipeline

```commandline (example)
cd /gpfs/data/cbc/sanghyun_lee/
nextflow_start
module load java/jdk-11

nextflow run /gpfs/data/cbc/pcao5/workflow_prototypes_on_OSCAR/workflows/rnaseq_homer.nf \
--samplesheet /gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F001_analysis/samplesheet_homer_sanghyun.csv \
--out_dir homer_out --reference_genome /gpfs/data/cbc/sanghyun_lee/mm10/star_idx4 \
--homer_data /gpfs/data/cbc/sanghyun_lee/homer/data --homer_config /gpfs/data/cbc/sanghyun_lee/homer/config.txt \
--erv_gtf /gpfs/data/cbc/sanghyun_lee/Mmus38.geve.m_v2_SRK.gtf
```

```commandline
cd $WORK_DIR
nextflow_start
module load java/jdk-11

nextflow run /gpfs/data/cbc/pcao5/workflow_prototypes_on_OSCAR/workflows/rnaseq_homer.nf \
--samplesheet $SAMPLESHEET \
--out_dir $OUTPUT_DIRECTORY --reference_genome $STAR_INDEX_DIRECTORY \
--homer_data $HOMER_DATA_DIRECTORY --homer_config $HOMER_CONFIG_FILE \
--erv_gtf $ERV_GTF
```

#### 4a. HOMER Samplesheet

```commandline
sample_id,read1,read2,group
S31,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F001/01.RawData/S31/S31_1.fq.gz,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F0
01/01.RawData/S31/S31_2.fq.gz,KO
S32,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F001/01.RawData/S32/S32_1.fq.gz,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F0
01/01.RawData/S32/S32_2.fq.gz,KO
S33,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F001/01.RawData/S33/S33_1.fq.gz,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F0
01/01.RawData/S33/S33_2.fq.gz,KO
S34,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F001/01.RawData/S34/S34_1.fq.gz,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F0
01/01.RawData/S34/S34_2.fq.gz,WT
S35,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F001/01.RawData/S35/S35_1.fq.gz,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F0
01/01.RawData/S35/S35_2.fq.gz,WT
S36,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F001/01.RawData/S36/S36_1.fq.gz,/gpfs/data/cbc/sanghyun_lee/X202SC22122741-Z01-F0
01/01.RawData/S36/S36_2.fq.gz,WT
```


### Appendix

### 2. Build the Docker
```bash
cd metadata
docker build -f RNASeq_Dockerfile -t rnaseq_pipeline --no-cache --platform linux/amd64 .
```

### 3. Test Run
```commandline
nextflow run workflows/rnaseq.nf --samplesheet metadata/samplesheet_rnaseq.csv --out_dir out
```