#! /usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.read && !params.sequencing_summary && !params.fast5) {
  error "Error: at least one input format (--read, --sequencing_summary, --fast5) must be enabled."
}

process guppy {
  container 'cowmoo/pycoqc:latest'

  input:
    file fast5s

  output:
    path "guppy_out/*.fastq", emit: reads
    path 'guppy_out/sequencing_summary.txt', emit: sequencing_summary

  """
  /ont-guppy-cpu/bin/guppy_basecaller -i ${fast5s} -s guppy_out -c dna_r9.4.1_450bps_hac.cfg --num_callers 8 --cpu_threads_per_caller 1
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

process pycoqc {
  container 'cowmoo/pycoqc:latest'

  input:
    file reads_summary

  output:
    path "pycoQC_output.html"

  """
  pycoQC -f ${reads_summary} -o pycoQC_output.html
  """
}

process kraken {
  container 'cowmoo/pycoqc:latest'
  containerOptions '-v /Users/test_user/minikraken2_v2_8GB_201904_UPDATE:/db'

  input:
    file reads

  output:
    path "out.txt"

  """
  /kraken2/kraken2 -db /db ${reads} --gzip-compressed --output out.txt --report report.txt
  """
}

process nanoPlot {
  container 'cowmoo/pycoqc:latest'

  input:
    file input

  output:
    path "summary-plots-log-transformed/*"

  script:
  if (params.reads_summary) {
   """
   NanoPlot --summary ${input} --loglength -o summary-plots-log-transformed
   """
  }
  else if (params.read) {
   """
   NanoPlot --fastq ${input} --loglength -o summary-plots-log-transformed
   """
  }
}

workflow {

  /* run guppy if fast5 directory is supplied */
  if (params.fast5s) {
    fast5 = file(params.fast5)
    guppy(fast5)
  }

  /* run taxonomical identification, either by guppy's called reads or user-suplied reads */
  read = null
  if (params.fast5s || params.read) {
    read = (params.fast5s) ? concat_reads(guppy.out.reads.collect()) : file(params.read)
    kraken(read)
  }

  /* run reads QC, either by the raw reads or user-supplied reads_summary.txt */
  if (params.reads_summary) {
    reads_summary = file(reads_summary)

    nanoPlot(reads_summary)
    pycoQC(reads_summary)
  } else if (params.fast5s) {
    reads_summary = file(guppy.out.reads_summary.collect())

    nanoPlot(reads_summary)
    pycoQC(reads_summary)
  } else if (read) {
    nanoPlot(read)
  }

}
