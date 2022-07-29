import nextflow.splitter.CsvSplitter

nextflow.enable.dsl=2

def fetchRunAccessions( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> run_accessions = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}

process downloadFiles {
  input:
  val id
  output:
  tuple val(id), path("${id}**.fastq")
  """
  fasterq-dump --split-3 ${id}
  """
}

process bowtie2 {
  label 'align'
  input:
  tuple val(sample), path(readsFastq)

  output:
  tuple val(sample), path("numReads.txt"), path("alignments*.sam")

  script:
  if(params.libraryLayout == 'single')
      """
      grep -c '^@' ${readsFastq} > numReads.txt

      ${params.bowtie2Command} \
        -x ${params.refdb} \
        -U ${readsFastq} \
        -S alignmentsSingle.sam 
      """
  else if(params.libraryLayout == 'paired')
      """
      grep -c '^@' ${sample}_1.fastq > numReads.txt

      ${params.bowtie2Command} \
        -x ${params.refdb} \
        -1 ${sample}_1.fastq \
        -2 ${sample}_2.fastq \
        -S alignmentsPaired.sam
      """
}

process alignmentStats {
  publishDir "${params.resultDir}/alignmentStats"
  label 'stats'
  input:
  tuple val(sample), path(numReadsPath), path(alignmentsSam)

  output:
  tuple val(sample), path("${sample}.alignmentStats.txt")

  script:
  """
  ${params.alignmentStatsCommand} ${alignmentsSam} > ${sample}.alignmentStats.txt
  """
}

process summarizeAlignments{
  publishDir "${params.resultDir}/summarizedAlignments"
  label 'postAlign'
  input:
  tuple val(sample), path(numReadsPath), path(alignmentsSam)

  output:
  path("${sample}.taxa.tsv")

  script:
  """
  ${params.summarizeAlignmentsCommand} \
    --input ${alignmentsSam} \
    --refdb-marker-to-taxon-path ${params.markerToTaxonPath} \
    --refdb-format eukprot \
    --output-type taxon_all \
    --num-reads \$(cat ${numReadsPath}) \
    --output ${sample}.taxa.tsv 
  """
}

process makeTsv {
  publishDir params.resultDir, mode: 'move', overwrite: true  
  label 'postAlign'

  input:
  file("*.taxa.tsv")

  output:
  file("${params.summaryColumn}.${params.summaryFormat}.tsv")

  script:
  """
  makeTsv.pl . .taxa.tsv ${params.summaryColumn} ${params.summaryFormat} > ${params.summaryColumn}.${params.summaryFormat}.tsv
  """
}

def postAlign(sample_numReadsPath_alignmentsSam) {
  alignmentStats(sample_numReadsPath_alignmentsSam)
  return summarizeAlignments(sample_numReadsPath_alignmentsSam)
}

workflow {
  if (params.downloadMethod == 'sra') {
    accessions = fetchRunAccessions(params.inputPath)
    ids = Channel.fromList(accessions)
    sample_reads = downloadFiles(ids)
  } else if (params.downloadMethod == 'local') {
    sample_reads = Channel.fromPath(params.inputPath).splitCsv(sep: "\t")
  }
  sample_numReads_alignments = bowtie2(sample_reads)
  xs = postAlign(sample_numReads_alignments)
  makeTsv(xs.collect())
}
