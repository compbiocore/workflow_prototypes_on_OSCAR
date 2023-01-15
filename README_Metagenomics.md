### Run Metagenomics Workflow

```bash
cd /gpfs/data/cbc/pcao5
nextflow -bg run workflows/metagenomics.nf --samplesheet metadata/samplesheet_metagenomics.csv --metadata metadata/fastqscreen.conf --outdir out
```