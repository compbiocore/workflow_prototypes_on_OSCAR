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
   """

   if (read2.exists)
     """
     fastqc ${read2}
     """
}

process trimmomatic {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.outdir"

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  input:
    tuple val(sample_id), file(read1), file(read2)

  output:
    path "*"
    file("fastq/${sample_id}_tr.fq.gz"), emit: fastq
    tuple val(sample_id), file(''), file("fastq/${sample_id}_tr.fq.gz"), file(null), emit: out

  script:
    """
     mkdir fastq logs
     TrimmomaticPE -threads 8 -trimlog logs/${sample_id}_trimmomatic_PE.log ${read1} ${read2} -baseout fastq/${sample_id}_tr.fq.gz ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:5:6:true SLIDINGWINDOW:10:25 MINLEN:50
    """
}

workflow PROCESS_SAMPLE {
    take:
        input_ch
    main:
        fastqc(input_ch)
        fastqc(trimmomatic(input_ch).out)
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
