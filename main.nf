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

process prepSingleSra {

  label 'prep'

  input:
  tuple val(sample), path(runAccession)

  output:
  tuple val(sample), path("${sample}.fastq")

  """
  gzip -d --force $runAccession
  """
}

process prepPairedSra {

  label 'prep'

  input:
  tuple val(sample), path(runAccession)

  output:
  tuple val(sample), path("${sample}_1.fastq"), path("${sample}_2.fastq")

  """
  gzip -d --force ${runAccession[0]} 
  gzip -d --force ${runAccession[1]}
  """
}

process bowtie2Single {
  label 'align'
  input:
  tuple val(sample), path(readsFastq)

  output:
  tuple val(sample), path("numReads.txt"), path("alignmentsSingle.sam")

  script:
  """
  grep -c '^@' ${readsFastq} > numReads.txt

  ${params.bowtie2Command} \
    -x ${params.refdb} \
    -U ${readsFastq} \
    -S alignmentsSingle.sam 
  """

}

process bowtie2Paired {
  label 'align'
  input:
  tuple val(sample), path(readsFastqR1), path(readsFastqR2)

  output:
  tuple val(sample), path("numReads.txt"), path("alignmentsPaired.sam")

  script:
  """
  grep -c '^@' ${readsFastqR1} > numReads.txt

  ${params.bowtie2Command} \
    -x ${params.refdb} \
    -1 ${readsFastqR1} \
    -2 ${readsFastqR2} \
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
    input = Channel.fromSRA( accessions, apiKey: params.apiKey, protocol: "http" )
  } else if (params.downloadMethod == 'local') {
    sample_reads = Channel.fromPath(params.inputPath).splitCsv(sep: "\t")
  }
  if(params.downloadMethod == 'sra' && params.libraryLayout == 'single'){
    sample_reads = prepSingleSra(input)
  } else if(params.downloadMethod == 'sra' && params.libraryLayout == 'paired') {
    sample_reads = prepPairedSra(input)
  }
  if(params.libraryLayout == 'single') {
    sample_numReads_alignments = bowtie2Single(sample_reads)
  } else if (params.libraryLayout == 'paired') {
    sample_numReads_alignments = bowtie2Paired(sample_reads)
  }
  xs = postAlign(sample_numReads_alignments)
  makeTsv(xs.collect())
}
