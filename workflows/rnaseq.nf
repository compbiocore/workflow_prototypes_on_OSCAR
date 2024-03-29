#! /usr/bin/env nextflow
nextflow.enable.dsl=2

params.reference_genome = "/gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/references/Oryctolagus_cuniculus.OryCun2.0_star_idx"
params.gtf = "/gpfs/data/cbc/koren_lab/elif_sengun_rnaseq_ffs/references/Oryctolagus_cuniculus.OryCun2.0.108.gtf"
params.htseq_multisample = false
params.sjdbGTFfile = 99
params.reference_genome_fasta = ""
params.fastqscreen_conf = "/gpfs/data/cbc/pcao5/workflow_prototypes_on_OSCAR/metadata/fastqscreen_mouse_rnaseq.conf"
params.qc_only = false
params.kraken = false

if (!params.samplesheet || !params.out_dir) {
  error "Error: Missing the samplesheet (--samplesheet) or output directory (--out_dir)."
}

process fastq_screen {

  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  time '12.h'

  cpus 8

  memory '100.GB'

  publishDir "$params.out_dir/qc/", mode: 'copy', overwrite: false

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  input:
    tuple val(sample_id), path(reads)

  output:
    path "${sample_id}_fastq_screen/*"

  """
  /FastQ-Screen-0.15.2/fastq_screen --aligner bwa --conf ${file(params.fastqscreen_conf)} --outdir ${sample_id}_fastq_screen ${reads}
  """
}

process multiqc {
  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  publishDir "$params.out_dir/qc", mode: 'copy', overwrite: false

  input:
   path(fastqcs)
   path(fastq_screens)

  output:
   path("*")

  script:
   """
    multiqc *_fastqc.zip *_screen.txt
   """
}

process multiqc_full {
  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  publishDir "$params.out_dir/qc", mode: 'copy', overwrite: false

  input:
   path(fastqcs)
   path(fastq_screens)
   path(qualimap)
   path(htseq_count)

  output:
   path("*")

  script:
   """
    multiqc *_fastqc.zip *_screen.txt */rnaseq_qc_results.txt *_htseq_counts
   """
}

process multiqc_erv {
  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  publishDir "$params.out_dir/qc", mode: 'copy', overwrite: false

  input:
   path(fastqcs)
   path(fastq_screens)
   path(qualimap)
   path(htseq_count)

  output:
   path("*")

  script:
   """
    multiqc *_fastqc.zip *_screen.txt */genome_results.txt *_htseq_counts
   """
}

process build_star_index {
  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  cpus 6

  memory '25.GB'

  time '6.h'

  output:
   path ("genome_idx")

  script:
   """
    STAR --runThreadN 6 \
         --runMode genomeGenerate \
         --genomeDir genome_idx \
         --genomeSAindexNbases 12 \
         --genomeFastaFiles ${file(params.reference_genome_fasta)} \
         --sjdbGTFfile ${file(params.gtf)} \
         --sjdbOverhang ${params.sjdbGTFfile}
   """
}



process qualimap {
  errorStrategy 'ignore'

  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  publishDir "$params.out_dir/qc/${sample_id}", mode: 'copy', overwrite: false

  memory '25.GB'

  time '6.h'

  input:
   tuple val(sample_id), file(bam), file(bam_index)

  output:
   path "${sample_id}"

  script:
   """
    qualimap rnaseq -gtf ${params.gtf} -bam ${bam} -outdir ${sample_id} --java-mem-size=25G
   """
}

process qualimap_erv {
  errorStrategy 'ignore'

  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  publishDir "$params.out_dir/qc/${sample_id}", mode: 'copy', overwrite: false

  memory '25.GB'

  time '6.h'

  cpus 5 

  input:
   tuple val(sample_id), file(bam), file(bam_index)

  output:
   path "${sample_id}"

  script:
   """
    qualimap bamqc -gff ${params.gtf} -bam ${bam} -outdir ${sample_id} --java-mem-size=25G
   """
}


process mark_duplicate {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir/", pattern: "*.dup.srtd.bam*", mode: 'copy', overwrite: false

  input:
    tuple val(sample_id), file(alignment)

  output:
    tuple val(sample_id), file("*.dup.srtd.bam"), file("*.dup.srtd.bam.bai"), emit: marked
    tuple path("*.dup.srtd.bam"), emit: bams

  script:
   """
     samtools index -b ${alignment}
     bammarkduplicates2 I=${alignment} O=${sample_id}.dup.srtd.bam M=${sample_id}.dup.metrics.txt
     samtools index -b ${sample_id}.dup.srtd.bam
   """
}

process fastqc {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir/qc", mode: 'copy', overwrite: false

  input:
    tuple val(sample_id), file(read1), file(read2)

  output:
    path "*_fastqc.*"

  script:
   if (read2.size() > 0)
     """
      fastqc ${read1}
      fastqc ${read2}
     """
   else
     """
      fastqc ${read1}
     """
}

process fastqc2 {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir/qc", mode: 'copy', overwrite: false

  input:
    tuple val(sample_id), path(reads)

  output:
    path "*_fastqc.*"

  script:
    if (reads.size() > 1)
     """
      fastqc ${reads[0]}
      fastqc ${reads[1]}
     """
   else
     """
      fastqc ${reads[0]}
     """
}

process kraken {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir", pattern: "*.txt", mode: 'copy', overwrite: false

  containerOptions '--bind /gpfs/data/cbc/pcao5/minikraken2_v2_8GB_201904_UPDATE:/db'

  memory '30.GB'

  time '2.h'

  memory '10.GB'

  input:
    tuple val(sample_id), file(read1), file(read2)

  output:
    path "*.txt"
  
  script:
   if (read2.size() > 0)
    """
     /kraken2-2.1.2/kraken2/kraken2 --paired -db /db ${read1} ${read2} --gzip-compressed --output ${sample_id}_kraken_log.txt --report ${sample_id}_kraken.txt
    """
   else
    """
     /kraken2-2.1.2/kraken2/kraken2 -db /db ${read1} --gzip-compressed --output ${sample_id}_kraken_log.txt --report ${sample_id}_kraken.txt
    """
}

process trimmomatic {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir", pattern: "*.fq.gz", mode: 'copy', overwrite: false

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  time '6.h'

  cpus 8

  memory '25.GB'

  tag "${sample_id}.trimmomatic"

  input:
    tuple val(sample_id), file(read1), file(read2)

  output:
    tuple val(sample_id), path("fastq/*P.fq.gz")

  script:
    if (read2.size() > 0)
     """
      mkdir fastq logs
      TrimmomaticPE -threads 8 -trimlog logs/${sample_id}_trimmomatic_PE.log ${read1} ${read2} -baseout fastq/${sample_id}_tr.fq.gz ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:5:6:true SLIDINGWINDOW:10:25 MINLEN:50
     """
    else
     """
      mkdir fastq logs
      TrimmomaticSE -threads 8 -trimlog logs/${sample_id}_trimmomatic_SE.log ${read1} fastq/${sample_id}_trP.fq.gz ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:5:6:true SLIDINGWINDOW:10:25 MINLEN:50
     """
}

process star {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir", pattern: "*.bam", mode: 'copy', overwrite: false

  time '6.h'

  cpus 16

  memory '75.GB'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs/:/gpfs/data/shared/databases/refchef_refs/'

  input:
    tuple val(sample_id), path(reads)
    path(reference_genome)

  output:
    tuple val(sample_id), file("*.sortedByCoord.out.bam")

  script:
    """
     STAR --genomeLoad NoSharedMemory --runThreadN 16 --outBAMsortingThreadN 12 --genomeDir ${reference_genome} \
          --quantMode GeneCounts --twopassMode Basic --sjdbGTFfile ${params.gtf} -outReadsUnmapped Fastx \
          --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --readFilesIn ${reads} \
          --outFileNamePrefix ${sample_id}
    """
}

process star_erv {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir", pattern: "*.bam", mode: 'copy', overwrite: false

  time '6.h'

  cpus 16

  memory '75.GB'

  containerOptions "--bind /gpfs/data/cbc:/gpfs/data/cbc"

  input:
    tuple val(sample_id), path(reads)
    path(reference_genome)

  output:
    tuple val(sample_id), file("*.sortedByCoord.out.bam")

  script:
    """
     STAR --runMode alignReads --runThreadN 16 --outBAMsortingThreadN 12 --genomeDir ${reference_genome} \
     --outFilterMultimapNmax 1000 --outFilterMismatchNmax 6 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
     --outFilterScoreMin 50 --readFilesIn ${reads} --readFilesCommand zcat --outFileNamePrefix ${sample_id} \
     --limitOutSAMoneReadBytes 200000 --outSAMtype BAM SortedByCoordinate
    """
}

process htseq_count_multisample {

  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir", mode: 'copy', overwrite: false

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  cpus 16

  input:
    path(bams)

  output:
   path "*"

  script:
   $/
    samtools index -M ${bams}
    htseq-count -s no -t exon -f bam -a 0 -r pos --additional-attr=gene_name --nonunique=all -i gene_id \
    --secondary-alignments=score ${bams} ${params.gtf} > htseq_counts
    echo ${bams} | sed -e '1s/^/gene gene_name /;s/\.dup.srtd.bam//g' |  tr ' ' \\t | cat - htseq_counts > tmpfile && mv tmpfile htseq_counts
   /$
}

process htseq_count {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir", mode: 'copy', overwrite: false

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  cpus 16

  time '5.h'

  memory '30.GB' 

  input:
   tuple val(sample_id), file(bam), file(bam_index)

  output:
   path "*_htseq_counts"

  script:
   """
    htseq-count -s no -t exon -f bam -a 0 -r pos --additional-attr=gene_name --nonunique=all -i gene_id \
    --secondary-alignments=score ${bam} ${params.gtf} > ${sample_id}_htseq_counts
   """
}

process feature_count {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir/expressions/", mode: 'copy', overwrite: false

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  memory '8.GB'

  time '5.h'

  input:
   tuple val(sample_id), file(bam), file(bam_index)

  output:
   path "*"

  script:
   """
    featureCounts -p -s 0 -M --fracOverlap 0.8 -O -a ${params.gtf} -o ${sample_id}.featureCounts.txt ${bam}
   """
}

process feature_count_erv {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir/expressions/", mode: 'copy', overwrite: false

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  memory '8.GB'

  time '5.h'

  input:
   tuple val(sample_id), file(bam), file(bam_index)

  output:
   path "*"

  script:
   """
    featureCounts -p -M --primary -s 1 -a ${params.gtf} -o ${sample_id}.featureCounts.txt ${bam}
   """
}

workflow PROCESS_SAMPLE {
    take:
        input_ch
        reference_genome

    main:
        fastqc(input_ch)

        if (params.kraken) {
            kraken(input_ch)
        }

        trimmed_reads = trimmomatic(input_ch)
        fastqcs = fastqc2(trimmed_reads).collect()
        fastqc_screens = fastq_screen(trimmed_reads).collect()

        if (params.qc_only) {
            multiqc(fastqcs, fastqc_screens)
        }

        mark_duplicate_bams = null

        if (!params.qc_only) {

            if (!params.erv) {
                marked_duplicates_bams = mark_duplicate(star(trimmed_reads, reference_genome))
                qualimaps = qualimap(marked_duplicates_bams.marked).collect()
                feature_count(marked_duplicates_bams.marked)
            } else {
                marked_duplicates_bams = mark_duplicate(star_erv(trimmed_reads, reference_genome))
                qualimaps = qualimap_erv(marked_duplicates_bams.marked).collect()
                feature_count_erv(marked_duplicates_bams.marked)
            }

            if (!params.htseq_multisample) {
                htseq_counts = htseq_count(marked_duplicates_bams.marked)
            }


            mark_duplicate_bams = marked_duplicates_bams.bams.collect()
            
            if (!params.erv) {
            	multiqc_full(fastqcs, fastqc_screens, qualimaps, htseq_counts)
            } else {
		multiqc_erv(fastqcs, fastqc_screens, qualimaps, htseq_counts)
	    }	
        }

    emit:
        mark_duplicate_bams
}

// Function to resolve files
def get_sample_info(LinkedHashMap sample) {
    read1  = sample.read1 ? file(sample.read1, checkIfExists: true) : null
    read2 = sample.read2 ? file(sample.read2, checkIfExists: true) : null

    return [ sample.sample_id, read1, read2 ]
}

workflow {
     reference_genome = params.reference_genome

     if (params.reference_genome_fasta != "") {
        reference_genome = build_star_index()
     }

     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .map { get_sample_info(it) }.set { samples_ch }

     PROCESS_SAMPLE (samples_ch, reference_genome)

     if (params.htseq_multisample && !params.qc_only) {
        htseq_count_multisample(PROCESS_SAMPLE.out.collect())
     }
}

workflow.onComplete {
  out_dir = workflow.launchDir.resolve(params.out_dir).toString()

  log.info("Pipeline completed at: $workflow.complete")
  log.info("Execution status: ${ workflow.success ? 'OK' : 'failed' }")

  File file = new File(out_dir + "/status.txt")
  if (workflow.errorReport) {
   file.write "Pipeline completed at: $workflow.complete\n${workflow.errorReport}"
  } else {
   file.write("Pipeline completed at: $workflow.complete\nExecution status: ${ workflow.success ? 'OK' : 'failed' }\n")
  }
}
