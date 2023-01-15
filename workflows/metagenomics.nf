#! /usr/bin/env nextflow
nextflow.enable.dsl=2

params.emudb = '/gpfs/data/cbc/pcao5/silva_db'

if (!params.samplesheet || !params.outdir || !params.fastqscreen_conf) {
  error "Error: Missing the samplesheet (--samplesheet), fastqscreen conf file (--fastqscreen_conf) or output directory (--outdir)."
}

process guppy {
  container 'cowmoo/pycoqc:latest'

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5), path(fastq_directory)

  output:
    tuple val(sample_id), file('guppy_out/sequencing_summary.txt'), file("guppy_out/*.fastq"), file(fast5), file(fastq_directory)

  script:
    if (fast5)
      """
      /ont-guppy-cpu/bin/guppy_basecaller -i ${fast5} -s guppy_out -c dna_r9.4.1_450bps_hac.cfg --num_callers 8 --cpu_threads_per_caller 1
      """
}

process concat_reads {
  debug true 

  container 'cowmoo/pycoqc:latest'

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5), path(fastq_directory)

  output:
    tuple val(sample_id), file(summary), file("reads.fastq.gz"), file(fast5), path(fastq_directory)

  script:
    if (fastq_directory)
      """
      cat ${fastq_directory}/*.fastq.gz > reads.fastq.gz
      """
}

process pycoQC {
  debug true

  container 'cowmoo/pycoqc:latest'

  publishDir "$params.outdir/$sample_id"

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5), path(fastq_directory)

  output:
    path "pycoQC_output.html"

  script:
    if (summary && summary.size() > 0)
      """
      pycoQC -f ${summary} -o pycoQC_output.html
      """
}

process fastq_screen {
  debug true

  container 'cowmoo/pycoqc:latest'

  memory '75.GB'

  time '6.h'

  cpus 4

  publishDir "$params.outdir/$sample_id/"

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5), path(fastq_directory)

  output:
    path "fastq_screen/*"

  """
  /FastQ-Screen-0.15.2/fastq_screen --aligner bwa --conf ${file(params.fastqscreen_conf)} --outdir fastq_screen ${reads}
  """
}

process emu {
  debug true
  
  cpus 6

  time '1.d'

  container 'cowmoo/pycoqc:latest'

  label 'OSCAR_large_job'

  publishDir "$params.outdir/$sample_id/"

  containerOptions "--bind $params.emudb:/emu_db"

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5), path(fastq_directory)

  output:
    path "emu/*_rel-abundance.tsv"

  script:
    if (reads)
        """
        source /opt/conda/bin/activate py37
        export EMU_DATABASE_DIR=/emu_db
        emu abundance --threads 6 ${reads} --output-dir emu
        """
}

process nanoPlot {
    debug true

    container 'cowmoo/pycoqc:latest'

    publishDir "$params.outdir/$sample_id/"

    input:
        tuple val(sample_id), file(summary), file(reads), file(fast5), path(fastq_directory)

    output:
        path "nanoplot/*"

    script:
        if (summary && summary.size() > 0)
            """
            NanoPlot --summary ${summary} --loglength -o nanoplot
            """
        else if (reads)
            """
            NanoPlot --fastq ${reads} --loglength -o nanoplot
            """
}

workflow CONCAT_FASTQ {
    take:
        fastq_directory_ch

    main:
        fastq_directory_ch.map { sample -> tuple(sample.sample_id, null, null, null, file(sample.fastq_directory)) } | concat_reads

    emit:
        concat_reads.out
}

workflow BASE_CALL {
    take:
        fast5_ch

    main:
        fast5_ch.map { sample -> tuple(sample.sample_id, null, null, file(sample.fast5), null) } | guppy

    emit:
        guppy.out
}

workflow PROCESS_SAMPLE {
    take:
        input_ch
    main:
        nanoPlot(input_ch)
	    emu(input_ch)
	    fastq_screen(input_ch)
    emit:
        fastq_screen.out
}

// Function to resolve files
def get_sample_info(LinkedHashMap sample) {
    summary = sample.read_summary ? file(sample.read_summary, checkIfExists: true) : null
    reads = sample.reads ? file (sample.reads, checkIfExists: true) : null
    fast5 = sample.fast5 ? file(sample.fast5, checkIfExists: true) : null

    return [ sample.sample_id, summary, reads, fast5, sample.fastq_directory ]
}

workflow {
     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .filter(row -> row.fast5)
            .set { fast5_ch }
     BASE_CALL (fast5_ch).set { base_called_ch }

     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .filter(row -> row.fastq_directory)
            .set { fastq_directory_ch }
     CONCAT_FASTQ (fastq_directory_ch).set{ combined_fastq_ch }

     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .filter(row -> (row.summary || row.reads) && (!row.fast5 || !row.fastq_directory))
            .map { get_sample_info(it) }.set { samples_ch }

     PROCESS_SAMPLE (samples_ch.concat(base_called_ch).concat(combined_fastq_ch))

     emit:
        PROCESS_SAMPLE.out
}

