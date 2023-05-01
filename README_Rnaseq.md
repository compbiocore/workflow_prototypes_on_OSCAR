### 1. Pull Down the Repo

If you haven't already pulled down the repo

```bash 
cd $HOME
git clone https://github.com/compbiocore/workflows_on_OSCAR.git
```

### 2. Modify the CSV File

Edit your CSV file e.g., `samplesheet_rnaseq.csv`

```commandline
## Fix the CSV file
sample_id,read1,read2
aged_1,/gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/00_fastq/B440-RA_R1_001.fastq.gz,/gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/00_fastq/B440-RA_R2_001.fastq.gz
young_1,/gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/00_fastq/F67-RA_R1_001.fastq.gz,/gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/00_fastq/F67-RA_R2_001.fastq.gz
```

### 3. Run the RNASeq

You can use `-bg` flag to run your workflow in the background (equivalent to `nohup` so you log off your interactive session).

```bash
cd $HOME/workflow_prototypes_on_OSCAR
nextflow_start

nextflow -bg run workflows/rnaseq.nf --samplesheet samplesheet_rnaseq.csv --out_dir rnaseq_test_out \
--reference_genome /gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/references/Oryctolagus_cuniculus.OryCun2.0_star_idx --gtf /gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/references/Oryctolagus_cuniculus.OryCun2.0.108.gtf
```

### 4. Check Your Output Folder

The outputs for samples should be organized by `expressions/$SAMPLE_NAME*`, `qc/$SAMPLE_NAME*`, `$SAMPLE_NAME*.bam`, `$SAMPLE_NAME_htseq_counts`.

```bash
cd $HOME/workflow_prototypes_on_OSCAR
ls -R rnaseq_test_out

rnaseq_test_out:
aged_1.dup.srtd.bam                  aged_1_htseq_counts  report.html           young_1.dup.srtd.bam.bai
aged_1.dup.srtd.bam.bai              expressions          timeline.html         young_1Aligned.sortedByCoord.out.bam
aged_1Aligned.sortedByCoord.out.bam  qc                   young_1.dup.srtd.bam  young_1_htseq_counts

rnaseq_test_out/expressions:
aged_1.featureCounts.txt  aged_1.featureCounts.txt.summary  young_1.featureCounts.txt  young_1.featureCounts.txt.summary

rnaseq_test_out/qc:
B440-RA_R1_001_fastqc.html  F67-RA_R1_001_fastqc.zip   aged_1_tr_1P_fastqc.zip    young_1_tr_1P_fastqc.zip
B440-RA_R1_001_fastqc.zip   F67-RA_R2_001_fastqc.html  aged_1_tr_2P_fastqc.html   young_1_tr_2P_fastqc.html
B440-RA_R2_001_fastqc.html  F67-RA_R2_001_fastqc.zip   aged_1_tr_2P_fastqc.zip    young_1_tr_2P_fastqc.zip
B440-RA_R2_001_fastqc.zip   aged_1                     young_1
F67-RA_R1_001_fastqc.html   aged_1_tr_1P_fastqc.html   young_1_tr_1P_fastqc.html

rnaseq_test_out/qc/aged_1:
aged_1

rnaseq_test_out/qc/aged_1/aged_1:
css  images_qualimapReport  qualimapReport.html  raw_data_qualimapReport  rnaseq_qc_results.txt

rnaseq_test_out/qc/aged_1/aged_1/css:
agogo.css        bgtop.png           doctools.js       jquery.js     qualimap_logo_small.png  up-pressed.png
ajax-loader.gif  comment-bright.png  down-pressed.png  minus.png     report.css               up.png
basic.css        comment-close.png   down.png          plus.png      searchtools.js           websupport.js
bgfooter.png     comment.png         file.png          pygments.css  underscore.js

rnaseq_test_out/qc/aged_1/aged_1/images_qualimapReport:
Coverage Profile Along Genes (High).png  Coverage Profile Along Genes (Total).png  Reads Genomic Origin.png
Coverage Profile Along Genes (Low).png   Junction Analysis.png                     Transcript coverage histogram.png

rnaseq_test_out/qc/aged_1/aged_1/raw_data_qualimapReport:
coverage_profile_along_genes_(high).txt  coverage_profile_along_genes_(low).txt  coverage_profile_along_genes_(total).txt

rnaseq_test_out/qc/young_1:
young_1

rnaseq_test_out/qc/young_1/young_1:
css  images_qualimapReport  qualimapReport.html  raw_data_qualimapReport  rnaseq_qc_results.txt

rnaseq_test_out/qc/young_1/young_1/css:
agogo.css        bgtop.png           doctools.js       jquery.js     qualimap_logo_small.png  up-pressed.png
ajax-loader.gif  comment-bright.png  down-pressed.png  minus.png     report.css               up.png
basic.css        comment-close.png   down.png          plus.png      searchtools.js           websupport.js
bgfooter.png     comment.png         file.png          pygments.css  underscore.js

rnaseq_test_out/qc/young_1/young_1/images_qualimapReport:
Coverage Profile Along Genes (High).png  Coverage Profile Along Genes (Total).png  Reads Genomic Origin.png
Coverage Profile Along Genes (Low).png   Junction Analysis.png                     Transcript coverage histogram.png

rnaseq_test_out/qc/young_1/young_1/raw_data_qualimapReport:
coverage_profile_along_genes_(high).txt  coverage_profile_along_genes_(low).txt  coverage_profile_along_genes_(total).txt
```

### 4. Logging and Debugging

The logs for samples will organized under `output_dir/logs/${TASK_NAME/`;
Each task will have its corresponding:
 - executed shell command `.command.sh.txt`
 - stdout `.command.err.log.txt`
 - stderr `.command.err.err.txt`

In addition the Slurm execution will also have its corresponding 
 - executed job submission command `.command.run.txt`
 - job sumbission stdout/stderr `.command.begin.txt`

Each task will also have its own symlink reference to the original Nextflow working directory.

e.g., 

```bash
young_1_trimmomatic.command.run.txt
young_1_trimmomatic.command.sh.txt
young_1_trimmomatic.command.begin.txt
young_1_trimmomatic.command.err.txt
young_1_trimmomatic.command.log.txt
young_1_trimmomatic.command.out.txt
young_1_trimmomatic.command.trace.txt
young_1_trimmomatic -> /gpfs/data/nextflow_working_dir/71/2ae13fa7baaf4276107c2b07f3a255 #nextflow working directory symlink

young_2_trimmomatic.command.run.text
young_2__trimmomatic.command.sh.txt
young_2_trimmomatic => /gpfs/data/nextflow_working_dir/75/2ae13ba7cdaf42r6101c3e07f4d289 #nextflow working directory symlink
...
```

Finally each workflow's output directory will have a `status.txt` that will be written whenever a workflow finishes to indicate if the workflow finished successfully or it has failed and what the error is. 

#### 4a. Workflow Success Message for `status.txt`

```bash
Pipeline completed at: 2023-04-30T18:47:51.829232-04:00
Execution status: OK
```

#### 4b. Example Workflow Error Message for `status.txt`

The log should pinpoint to the failed job that caused the entire workflow to fail, its stdout, stderr and exit code; and the path to its Nextflow working directory, e.g.:
```bash
ipeline completed at: 2023-04-30T22:42:32.144450-04:00
Error executing process > 'sayHello'

Caused by:
  Process `sayHello` terminated with an error exit status (1)

Command executed:

  echo 'Hello world!'
  exit 1

Command exit status:
  1

Command output:
  Hello world!

Command wrapper:
  ## SLURM PROLOG ###############################################################
  ##    Job ID : 9846742
  ##  Job Name : nf-sayHello
  ##  Nodelist : node1323
  ##      CPUs : 4
  ##  Mem/Node : 4096 MB
  ## Directory : /gpfs/nextflow_publishdir/work/1f/50b11428d9238b52adccd961788510
  ##   Job Started : Sun Apr 30 22:42:29 EDT 2023
  ###############################################################################
  Hello world!

Work dir:
  /gpfs/nextflow_publishdir/work/1f/50b11428d9238b52adccd961788510

Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`
```

### Appendix:

### A. Build the Docker
```bash
cd metadata
docker build -f RNASeq_Dockerfile -t rnaseq_pipeline --no-cache --platform linux/amd64 .
```

### B. Test Quickly Workflow
```bash
nextflow run workflows/rnaseq.nf --samplesheet metadata/samplesheet_rnaseq_short.csv --out_dir rnaseq_test_out --reference_genome /gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/references/Oryctolagus_cuniculus.OryCun2.0_star_idx --gtf /gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/references/Oryctolagus_cuniculus.OryCun2.0.108.gtf
```