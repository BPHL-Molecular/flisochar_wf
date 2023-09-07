#!/usr/bin/env nextflow

/*
========================================================================================
                                     Flisochar
========================================================================================
Florida Isolate Characterization Pipeline using hybrid-assembly. Started 2022-08-22.
 #### Homepage / Documentation
 https://github.com/
 #### Author
 Tassy J. Bazile <tassy.bazile@flhealth.gov> - https://github.com/
----------------------------------------------------------------------------------------
*/

// Help message

def helpMessage() {
    log.info"""
    >>
    *** Florida Isolate Characterization Pipeline ***

    *** Flisochar ***

    Business Logic:
    This pipeline takes short and long-read fastq samples, generates qc reads, hybrid assemblies, identifies species, annotates genomes, seeks antimicrobial resistance genes, calculates average nucleotide identity using ANI.

    ========================================================================================================================================================================================
    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run flisochar_wf.nf --lr_type <pcb or ont> --lreads '*.fastq.gz' --sreads '*_{1,2}.fastq.gz' --asb_tool --genomeSize --outdir 'output path'
    Mandatory arguments:
        --lreads        Path to long-read fastq files
        --sreads        Path to short-read fastq files
        --asb_tool      One of the following assemblers: [canu, dragonflye, unicycler]
        --genomeSize    genome size when required, used with canu(ex. 3.5m) or dragonflye(ex. 3.5M)
        --outdir        Chosen name of the output directory
    Options:
        --kradb         Kraken maxi database path (if you want to use it)
        --krdbName      Kraken maxi dabase name in the path
         --threads      Number of cores
         -resume        If a pipeline workflow has been interrupted or stopped, this allows the workflow to pick up from the point it was interrupted.
        --lr_type       Types of long read technology either (ont, pcb), default: pcb
        --email         Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        -name           Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    Under conditions:
        --genomeSize    genome size not needed if run unicycler assembler
        <<
        """.stripIndent()

}

//Show help message

if (params.help){
    helpMessage()
    exit 0
}

// Input parameters

nextflow.enable.dsl = 2

params.lr_type = ""
params.sreads = ""
params.lreads = ""
params.outdir = "$baseDir/flisochar_results"
params.asb_tool = "dragonflye"
params.help = false
params.email = false
params.name = false
params.threads = false
params.genomeSize = 0
params.kradb = false
params.krdbName = false

// Whether the user specifies run name.
// bonus effet of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Input parameters and Validation
if (!params.lreads | !params.sreads){
    log.info"""
###################################################################
Please specify --lreads and --sreads as appropriate.
###################################################################
"""
    helpMessage()
    exit 0
}

// Create channels for input short and long-read files (, checkIfExists: true an alternative to ifEmpy)

lreads_ch = Channel
        .fromPath(params.lreads)
        .ifEmpty {exit 1, "Cannot find any long reads matching: ${params.lreads}\nNB: Path needs to be enclosed in quotes!"}
        .map { file -> tuple(file.simpleName, file) }
threads = Channel.value( params.threads )

sreads_ch = Channel
        .fromFilePairs(params.sreads, checkIfExists: true)
meta_ch = Channel.fromPath("$baseDir/metad/*.yaml").buffer(size: 2) // metadata channel to run pgap

// prints to the screen and to the log

log.info "========================================="
log.info " Flisochar v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Short Reads']  = params.sreads
summary['Long Reads']   = params.lreads
summary['Output dir']   = params.outdir
summary['Assembly tool'] = params.asb_tool
summary['Genome Size'] = params.genomeSize
summary['Long-read Technology']   = params.lr_type
summary['Kraken Maxi Database'] = params.kradb
summary['Name of Kraken Maxi Database'] = params.krdbName
summary['Working dir']  = workflow.workDir
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
................................................................
                        CONFIGURATION
................................................................
*/

// Show help message when not correct input

if (params.help) exit 0, helpMessage()

// Define valid assembly sofware tools
validAssemblyTools = ['dragonflye', 'unicycler', 'canu']
//value for parameters
asb_tools = params.asb_tool.tokenize(',') // replacing by a random value

if (!asb_tools.every{validAssemblyTools.contains(it)}) {
    log.info "Wrong execution asb_tool, should be one of " + validAssemblyTools
    exit 1
}

// Importing Module for each technology (ONT and PacBio)
include { ontflow } from './flisochar_modules/flisochar_ont.nf'
include { pcbflow } from './flisochar_modules/flisochar_pb_v3.nf'


workflow{
    if (params.lr_type == 'ont') {
        ontflow()
    }
    if (params.lr_type == 'pcb') {
        pcbflow()
    }
}

/*
workflow{
    if (params.lr_type == 'ont') {
        ontflow()
    }
    else{ 
        pcbflow()
    }
}
*/

