#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process downloadFiles {
  container = 'veupathdb/bowtiemapping'
  input:
    val id

  output:
    tuple val(id), path("${id}**.fastq")

  script:
    template 'downloadFiles.bash'
}


process bowtie2 {
  label 'align'

  input:
    tuple val(sample), path(readsFastq)

  output:
    tuple val(sample), path("numReads.txt"), path("alignments*.sam")

  script:
    if(params.libraryLayout.toLowerCase() == 'single')
      template 'bowtieSingle.bash'
    else if(params.libraryLayout.toLowerCase() == 'paired')
      template 'bowtiePaired.bash'
}


process alignmentStats {
  publishDir "${params.resultDir}/alignmentStats"

  label 'stats'

  input:
    tuple val(sample), path(numReadsPath), path(alignmentsSam)

  output:
    tuple val(sample), path("${sample}.alignmentStats.txt")

  script:
    template 'alignmentStats.bash'
}


process summarizeAlignments{
  publishDir "${params.resultDir}/summarizedAlignments"

  label 'postAlign'

  input:
    tuple val(sample), path(numReadsPath), path(alignmentsSam)

  output:
    path("${sample}.taxa.tsv")

  script:
    template 'summarizeAlignments.bash'
}


process makeTsv {
  publishDir params.resultDir, mode: 'move', overwrite: true  

  label 'postAlign'

  input:
    file("*.taxa.tsv")

  output:
    file("${params.summaryColumn}.${params.summaryFormat}.tsv")

  script:
    template 'makeTsv.bash'
}


def postAlign(sample_numReadsPath_alignmentsSam) {

  alignmentStats(sample_numReadsPath_alignmentsSam)
  return summarizeAlignments(sample_numReadsPath_alignmentsSam)

}


workflow sra {
  take:
    accessions

  main:
    ids = Channel.fromList( accessions )
    sample_reads = downloadFiles( ids )
    sample_numReads_alignments = bowtie2( sample_reads )
    xs = postAlign( sample_numReads_alignments )
    makeTsv(xs.collect())
}


workflow local {
  take:
    sample_reads

  main:
    sample_numReads_alignments = bowtie2( sample_reads )
    xs = postAlign( sample_numReads_alignments )
    makeTsv(xs.collect())
}
