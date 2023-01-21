#! /usr/bin/env nextflow
nextflow.enable.dsl=2

params.reference_genome = "/gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/references/Oryctolagus_cuniculus.OryCun2.0_star_idx"
params.gtf = "/gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/references/Oryctolagus_cuniculus.OryCun2.0_star_idx"


if (!params.samplesheet || !params.out_dir) {
  error "Error: Missing the samplesheet (--samplesheet) or output directory (--out_dir)."
}

process fastqc {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir/qc"

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

process fastqc2 {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir/qc"

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

process trimmomatic {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir"

  time '6.h'

  cpus 8

  input:
    tuple val(sample_id), file(read1), file(read2)

  output:
    path "*"
    tuple val(sample_id), file("fastq/${sample_id}_tr_1U.fq.gz"), file("fastq/${sample_id}_tr_2U.fq.gz"), emit: fastq_out

  script:
    """
     mkdir fastq logs
     TrimmomaticPE -threads 8 -trimlog logs/${sample_id}_trimmomatic_PE.log ${read1} ${read2} -baseout fastq/${sample_id}_tr.fq.gz ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:5:6:true SLIDINGWINDOW:10:25 MINLEN:50
    """
}

process star {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir"

  time '6.h'

  cpus 16

  memory '75.GB'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  input:
    tuple val(sample_id), file(read1), file(read2)

  script:
    """
     STAR --genomeLoad NoSharedMemory --runThreadN 16 --outBAMsortingThreadN 12 --genomeDir ${params.reference_genome} \
          --quantMode GeneCounts --twopassMode Basic --sjdbGTFfile ${params.gtf} -outReadsUnmapped Fastx \
          --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --readFilesIn ${read1} ${read2} \
          --outFileNamePrefix ${sample_id}
    """
}

workflow PROCESS_SAMPLE {
    take:
        input_ch
    main:
        fastqc(input_ch)
        trimmed_fastqs = trimmomatic(input_ch).fastq_out.collect()
        fastqc2(trimmed_fastqs)
        star(trimmed_fastqs)
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
