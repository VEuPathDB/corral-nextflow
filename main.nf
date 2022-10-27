#!/usr/bin/env nextflow

nextflow.enable.dsl=2
import nextflow.splitter.CsvSplitter

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

//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(params.downloadMethod.toLowerCase() == 'sra') {
  if (params.inputPath) {
    accessions = fetchRunAccessions( params.inputPath )
  }
  else {
    throw new Exception("Missing params.fastaSubsetSize")
  }
}
else if (params.downloadMethod.toLowerCase() == 'local') {
  sample_reads = Channel.fromFilePairs( params.inputPath + "/*_{1,2}.fastq" )
}
else {
  throw new Exception("Invalid value for params.downloadMethod")
}

//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { sra } from './modules/corral.nf'
include { local } from './modules/corral.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  if (params.downloadMethod.toLowerCase() == 'sra') {
    sra(accessions)
  }
  else if (params.downloadMethod.toLowerCase() == 'local') {
    local(sample_reads)
  }
}


