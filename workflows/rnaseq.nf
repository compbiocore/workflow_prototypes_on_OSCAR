#! /usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.samplesheet || !params.out_dir) {
  error "Error: Missing the samplesheet (--samplesheet) or output directory (--out_dir)."
}

process fastqc {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.outdir/qc"

  input:
    tuple val(sample_id), file(read1), file(read2)

  output:
    path "out/*"

  script:
    """
     fastqc -o out ${fastq} ${read1}
     fastqc -o out ${fastq} ${read2}
    """
}

workflow PROCESS_SAMPLE {
    take:
        input_ch
    main:
        fastqc(input_ch)
    emit:
        fastqc.out
}

workflow {
     Channel.fromPath(params.samplesheet).splitCsv(header:true).set { samples_ch }

     PROCESS_SAMPLE (samples_ch)

     emit:
        PROCESS_SAMPLE.out
}
