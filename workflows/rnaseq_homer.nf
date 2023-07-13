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

  time '6.h'

  cpus 6

  memory '25.GB'

  publishDir "$params.out_dir/qc/", mode: 'copy', overwrite: false

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  input:
    tuple val(sample_id), path(reads)

  output:
    path "${sample_id}_fastq_screen/*"

  """
  /FastQ-Screen-0.15.2/fastq_screen --aligner bwa --conf ${file(params.fastqscreen_conf)} --outdir ${sample_id}_fastq_screen ${reads} --subset 0
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
  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc --bind /gpfs/data/shared/databases/refchef_refs:/gpfs/data/shared/databases/refchef_refs'

  publishDir "$params.out_dir/qc/${sample_id}", mode: 'copy', overwrite: false

  memory '25.GB'

  time '6.h'

  input:
   tuple val(sample_id), file(bam), file(bam_index)

  output:
   path "*"

  script:
   """
    qualimap rnaseq -gtf ${params.gtf} -bam ${bam} -outdir ${sample_id}
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
    path "*"

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
    path "*"

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

  input:
    tuple val(sample_id), file(read1), file(read2), val(group)

  output:
    tuple val(sample_id), path("fastq/*P.fq.gz"), val(group)

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

process star_strict {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir", pattern: "*.bam", mode: 'copy', overwrite: false

  time '6.h'

  cpus 16

  memory '75.GB'

  containerOptions "--bind /gpfs/data/cbc:/gpfs/data/cbc --bind $params.homer_config:/homer/config.txt"

  input:
    tuple val(sample_id), path(reads), val(group)
    path(reference_genome)

  output:
    tuple path("${sample_id}"), emit: tag_directory
    tuple val(group), emit: groups

  script:
    """
     STAR --runMode alignReads --runThreadN 8 --genomeDir ${reference_genome} --outFilterMultimapNmax 1 \
     --readFilesIn ${reads} --readFilesCommand zcat --outFileNamePrefix ${sample_id}.Strict.mapped_to_mm10

     /homer/bin/makeTagDirectory ${sample_id} ${sample_id}.Strict.mapped_to_mm10Aligned.out.sam
    """
}

process star_loose {
  container 'cowmoo/rnaseq_pipeline:latest'

  publishDir "$params.out_dir", pattern: "*.bam", mode: 'copy', overwrite: false

  time '6.h'

  cpus 16

  memory '75.GB'

  containerOptions "--bind /gpfs/data/cbc:/gpfs/data/cbc --bind $params.homer_config:/homer/config.txt"

  input:
    tuple val(sample_id), path(reads), val(group)
    path(reference_genome)

  output:
    tuple path("${sample_id}"), emit: tag_directory
    tuple val(group), emit: groups

  script:
    """
     STAR --runMode alignReads --runThreadN 8 --genomeDir ${reference_genome} --outFilterMultimapNmax 1 \
     --readFilesIn ${reads} --readFilesCommand zcat --outFileNamePrefix ${sample_id}.Loose.mapped_to_mm10

     /homer/bin/makeTagDirectory ${sample_id} ${sample_id}.Loose.mapped_to_mm10Aligned.out.sam
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
   path "*"

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

process analyze_erv_repeats {
  container 'cowmoo/rnaseq_pipeline:latest'

  containerOptions "--bind $params.homer_data:/homer/data/ --bind $params.homer_config:/homer/config.txt"

  publishDir "$params.out_dir/erv/", mode: 'copy', overwrite: false

  memory '30.GB'

  time '5.h'

  input:
   path(tag_directory_path)
   val(groups)
   path(erv_gtf)

  output:
   file("countTable.Miner_Locus.txt")
   file("getDiff.Miner_Locus.txt")

  script:
   """
    export PATH=/homer/bin/:$PATH
    perl /homer/bin/analyzeRepeats.pl ${erv_gtf} mm10 -count genes -noadj -d ${tag_directory_path} > countTable.Miner_Locus.txt
    perl /homer/bin/getDiffExpression.pl countTable.Miner_Locus.txt ${groups.join(' ')} > getDiff.Miner_Locus.txt
   """
}

workflow PROCESS_SAMPLE {
    take:
        input_ch
        reference_genome

    main:
        trimmed_reads = trimmomatic(input_ch)

        strict_tag_directories = null
        group_labels = null

        if (!params.qc_only) {
            strict_bams = star_strict(trimmed_reads, reference_genome)
            loose_bams = star_loose(trimmed_reads, reference_genome)

            strict_tag_directories = strict_bams.tag_directory.collect()
            group_labels = strict_bams.groups.collect()
        }

    emit:
	    strict_tags = strict_tag_directories
	    groups = group_labels
}

// Function to resolve files
def get_sample_info(LinkedHashMap sample) {
    read1  = sample.read1 ? file(sample.read1, checkIfExists: true) : null
    read2 = sample.read2 ? file(sample.read2, checkIfExists: true) : null

    return [ sample.sample_id, read1, read2, sample.group ]
}

workflow {
     reference_genome = params.reference_genome

     if (params.reference_genome_fasta != "") {
        reference_genome = build_star_index()
     }

     Channel.fromPath(params.samplesheet).splitCsv(header:true)
            .map { get_sample_info(it) }.set { samples_ch }

     PROCESS_SAMPLE (samples_ch, reference_genome)

     analyze_erv_repeats(PROCESS_SAMPLE.out.strict_tags.collect(), PROCESS_SAMPLE.out.groups.collect(), params.erv_gtf)

}
