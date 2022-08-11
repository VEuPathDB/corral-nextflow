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

if(params.downloadMethod == 'sra') {
  if (params.inputPath) {
    accessions = fetchRunAccessions( params.inputPath )
  }
  else {
    throw new Exception("Missing params.fastaSubsetSize")
  }
}
else if (params.downloadMethod == 'local') {
  if (params.localFileLocation) {
  sample_reads = Channel.fromFilePairs( params.localFileLocation + "/*_{1,2}.fastq" )
  }
  else {
    throw new Exception("Missing params.localFileLocation")
  }
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
  if (params.downloadMethod == 'sra') {
    sra(accessions)
  }
  else if (params.downloadMethod == 'local') {
    local(sample_reads)
  }
}


