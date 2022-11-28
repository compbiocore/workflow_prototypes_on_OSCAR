#! /usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.samplesheet) {
  error "Error: Missing the samplesheet."
}

process guppy {
  container 'cowmoo/pycoqc:latest'

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5)

  output:
    tuple val(sample_id), file('guppy_out/sequencing_summary.txt'), file("guppy_out/*.fastq"), file(fast5)

  script:
    if (fast5)
      """
      /ont-guppy-cpu/bin/guppy_basecaller -i ${fast5} -s guppy_out -c dna_r9.4.1_450bps_hac.cfg --num_callers 8 --cpu_threads_per_caller 1
      """
}

process concat_reads {
  container 'cowmoo/pycoqc:latest'

  input:
    path reads

  output:
    file "reads.fastq", emit: read

  """
  cat ${reads} > reads.fastq
  """
}

process pycoQC {
  debug true

  container 'cowmoo/pycoqc:latest'

  publishDir "$params.outdir/$sample_id"

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5)

  output:
    path "pycoQC_output.html"

  script:
    if (summary)
      """
      pycoQC -f ${summary} -o pycoQC_output.html
      """
}

process kraken {
  debug true

  label 'OSCAR_large_job'

  container 'cowmoo/pycoqc:latest'

  publishDir "$params.outdir/$sample_id/kraken"

  containerOptions '--bind /users/test_user/nextflow_scratch/minikraken2_v2_8GB_201904_UPDATE:/db'

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5)

  output:
    path "out.txt"

  """
  /kraken2/kraken2 -db /db ${reads} --gzip-compressed --output out.txt --report report.txt
  """
}

process silva_minimap2 {
  debug true

  container 'cowmoo/pycoqc:latest'

  input:
    tuple val(sample_id), file(summary), file(reads), file(fast5)

  output:
    path "aln.sam"

  containerOptions '-v /Users/paulcao/ribogrove:/db'

  """
  /minimap2-2.24_x64-linux/minimap2 -t 5 -ax map-ont /db/ribogrove_8.214_sequences.fasta ${reads} > aln.sam
  """
}

process nanoPlot {
    debug true

    container 'cowmoo/pycoqc:latest'

    publishDir "$params.outdir/$sample_id"

    input:
        tuple val(sample_id), file(summary), file(reads), file(fast5)

    output:
        path "nanoplot/*"

    script:
        if (summary)
            """
            NanoPlot --summary ${summary} --loglength -o nanoplot
            """
        else if (reads)
            """
            NanoPlot --fastq ${reads} --loglength -o nanoplot
            """
}

/*
    organize each sample by their directory e.g.,
    WT/nanoPlot
    WT/pycoQC
    WT/kraken
*/

/*
process handoff {
}*/

workflow BASE_CALL {
    take:
        fast5_ch

    main:
        fast5_ch.map { sample -> tuple(sample.sample_id, null, null, file(sample.fast5)) } | guppy

    emit:
        guppy.out
}

workflow PROCESS_SAMPLE {
    take:
        input_ch
    main:
        nanoPlot(input_ch)
        pycoQC(input_ch)
        kraken(input_ch)
        /*silva_minimap2(input_ch)*/
    emit:
        nanoPlot.out
        /*kraken.out
        pycoQC.out
        /*silva_minimap2.out*/
}

// Function to resolve files
def get_sample_info(LinkedHashMap sample) {
    summary = sample.read_summary ? file(sample.read_summary, checkIfExists: true) : null
    reads = sample.reads ? file (sample.reads, checkIfExists: true) : null
    fast5 = sample.fast5 ? file(sample.fast5, checkIfExists: true) : null

    return [ sample.sample_id, summary, reads, fast5 ]
}

workflow {
     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .filter(row -> row.fast5)
            .set { fast5_ch }
     BASE_CALL (fast5_ch).set { base_called_ch }

     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .filter(row -> (row.summary || row.reads) && !row.fast5)
            .map { get_sample_info(it) }.set { samples_ch }

     PROCESS_SAMPLE (samples_ch.concat(base_called_ch))

     emit:
        PROCESS_SAMPLE.out
}
