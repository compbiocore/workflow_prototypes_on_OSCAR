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

  containerOptions '--bind /gpfs/data/cbc:/gpfs/data/cbc'

  publishDir "$params.out_dir/qc/${sample_id}", mode: 'copy', overwrite: false

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

process seqtk {
  container 'https://depot.galaxyproject.org/singularity/seqtk%3A1.3--h7132678_4'

  cpus 2

  memory '10.GB'

  input:
    tuple val(sample_id), file(read1), file(read2)

  output:
    tuple val(sample_id), path("sub_sampled.fasta")

  script:
   """
    seqtk sample -s 123 ${read1} 1000 > sub_sampled.fq.gz
    seqtk seq -a sub_sampled.fq.gz > sub_sampled.fasta
   """
}

process blastNR {
  container 'https://depot.galaxyproject.org/singularity/blast%3A2.13.0--hf3cf87c_0'

  containerOptions '-B /gpfs/data/shared/databases/refchef_refs/'

  cpus 6

  memory '25.GB'

  input:
    tuple val(sample_id), path(fasta)

  output:
    tuple val(sample_id), path("blast.xml")

  script:
   """
    blastn -num_threads 6 -query ${fasta} -db /gpfs/data/shared/databases/refchef_refs/nt_db/blast_db/nt -out blast.xml -outfmt 5
   """
}

process megan_process {
  container 'https://depot.galaxyproject.org/singularity/megan%3A6.24.20--h9ee0642_0'

  publishDir "$params.out_dir/", mode: 'copy', overwrite: false

  cpus 2

  memory '14.GB'

  input:
    tuple val(sample_id), path(xml)

  output:
    tuple val(sample_id), path("${sample_id}.out")

  script:
   """
    blast2rma -i ${xml} --format BlastXML --out out.rma
    rma2info -i out.rma -c2c Taxonomy -n -s -r -o ${sample_id}.out
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

workflow PROCESS_SAMPLE {
    take:
        input_ch

    main:
        kraken(input_ch)
        sampled_reads = seqtk(input_ch)
        blasts = blastNR(sampled_reads)
        megan_outs = megan_process(blasts)
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

}
