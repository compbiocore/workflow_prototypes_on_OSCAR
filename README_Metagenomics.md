### 1. Pull Down the Repo

If you haven't already pulled down the repo

```bash 
cd $HOME
git clone https://github.com/compbiocore/workflows_on_OSCAR.git
```

### 2. Download and Populate the Data
```bash 
module load sratoolkit
mkdir SRR19309654
cd SRR19309654/ && fastq-dump SRR19309654

mkdir SRR19309656
cd SRR19309656/ && fastq-dump SRR19309656
```

### 3. Modify the CSV File
Edit your CSV file
```commandline
## Fix the CSV file
sample_id,read_summary,reads,fast5,fastq_directory
SAMN28789430,,,,/gpfs/data/cbc/pcao5/workflow_prototypes_on_OSCAR/data/SRR19309654/
SAMN28789428,,,,/gpfs/data/cbc/pcao5/workflow_prototypes_on_OSCAR/data/SRR19309656/
```

### 4. Run Metagenomics Workflow

You can use `-bg` flag to run your workflow in the background (equivalent to `nohup` so you log off your interactive session).

```bash
cd $HOME/workflow_prototypes_on_OSCAR
nextflow_start
nextflow -bg run workflows/metagenomics.nf --samplesheet metadata/samplesheet_metagenomics.csv --fastqscreen_conf metadata/fastqscreen.conf --out_dir out --outdir out
```

### 5. Check Your Output Folder

The samples should be organized by `$SAMPLE_NAME/fastq_screen`, `$SAMPLE_NAME/nanoplot`, `$SAMPLE_NAME/emu`.

```bash
cd $HOME/workflow_prototypes_on_OSCAR
ls -R out

out:
SAMN28789428  SAMN28789430  report.html  timeline.html

out/SAMN28789428:
emu  fastq_screen  nanoplot

out/SAMN28789428/emu:
reads.fastq_rel-abundance.tsv

out/SAMN28789428/fastq_screen:
reads_screen.html  reads_screen.png  reads_screen.txt

out/SAMN28789428/nanoplot:
LengthvsQualityScatterPlot_dot.html            Non_weightedHistogramReadlength.html
LengthvsQualityScatterPlot_dot.png             Non_weightedHistogramReadlength.png
LengthvsQualityScatterPlot_kde.html            Non_weightedLogTransformed_HistogramReadlength.html
LengthvsQualityScatterPlot_kde.png             Non_weightedLogTransformed_HistogramReadlength.png
LengthvsQualityScatterPlot_loglength_dot.html  WeightedHistogramReadlength.html
LengthvsQualityScatterPlot_loglength_dot.png   WeightedHistogramReadlength.png
LengthvsQualityScatterPlot_loglength_kde.html  WeightedLogTransformed_HistogramReadlength.html
LengthvsQualityScatterPlot_loglength_kde.png   WeightedLogTransformed_HistogramReadlength.png
NanoPlot-report.html                           Yield_By_Length.html
NanoPlot_20230115_1844.log                     Yield_By_Length.png
NanoStats.txt

out/SAMN28789430:
emu  fastq_screen  nanoplot

out/SAMN28789430/emu:
reads.fastq_rel-abundance.tsv

out/SAMN28789430/fastq_screen:
reads_screen.html  reads_screen.png  reads_screen.txt

out/SAMN28789430/nanoplot:
LengthvsQualityScatterPlot_dot.html            Non_weightedHistogramReadlength.html
LengthvsQualityScatterPlot_dot.png             Non_weightedHistogramReadlength.png
LengthvsQualityScatterPlot_kde.html            Non_weightedLogTransformed_HistogramReadlength.html
LengthvsQualityScatterPlot_kde.png             Non_weightedLogTransformed_HistogramReadlength.png
LengthvsQualityScatterPlot_loglength_dot.html  WeightedHistogramReadlength.html
LengthvsQualityScatterPlot_loglength_dot.png   WeightedHistogramReadlength.png
LengthvsQualityScatterPlot_loglength_kde.html  WeightedLogTransformed_HistogramReadlength.html
LengthvsQualityScatterPlot_loglength_kde.png   WeightedLogTransformed_HistogramReadlength.png
NanoPlot-report.html                           Yield_By_Length.html
NanoPlot_20230115_1844.log                     Yield_By_Length.png
NanoStats.txt
```