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
    path "*"

  script:
    """
     fastqc ${read1}
     fastqc ${read2}
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

// Function to resolve files
def get_sample_info(LinkedHashMap sample) {
    read1  = sample.read1 ? file(sample.read1, checkIfExists: true) : null
    read2 = sample.read2 ? file(sample.read2, checkIfExists: true) : null

    return [ sample.sample_id, read1, read2 ]
}

workflow {
     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .map { get_sample_info(it) }.set { samples_ch }

     PROCESS_SAMPLE (samples_ch)

     emit:
        PROCESS_SAMPLE.out
}
