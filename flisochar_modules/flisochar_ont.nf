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
    
    ***	Flisochar ***
    
    Business Logic:
    This pipeline takes short and long-read fastq samples, generates qc reads, hybrid assemblies, identifies species, annotates genomes, seeks antimicrobial resistance genes, calculates average nucleotide identity using ANI.

    ========================================================================================================================================================================================
    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run flisochar_ont.nf --lreads '*.fastq.gz' --sreads '*_{1,2}.fastq.gz' --outdir 'output path'
    Mandatory arguments:
        --lreads        Path to long-read fastq files
        --sreads        Path to short-read fastq files
        --outdir        Chosen name of the output directory
    Options:
        --asb_tool      One of the following assemblers: [canu, dragonflye, unicycler]
        --threads       Number of cores
        --kradb         Kraken maxi database path (if you want to use it)
        --krdbName      Kraken maxi dabase name in the path
         -resume        If a pipeline workflow has been interrupted or stopped, this allows the workflow to pick up from the point it was interrupted.
        --email         Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        -name           Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    Under condition:
        --genomeSize    genome size when required, used with canu(ex. 3.5m) or dragonflye(ex. 3.5M)assembler 
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
params.sreads = ""
params.lreads = ""
params.outdir = "$baseDir/flisochar_results"
params.asb_tool = "dragonflye"
params.help = false
params.email = false
params.name = false
params.threads = (6)
params.genomeSize = "" // 09/06/23
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
        .map { file -> tuple(file.baseName.strip('.fastq'), file) }
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
summary['Assembly tool']   = params.asb_tool
summary['Genome Size']   = params.genomeSize 
summary['Working dir']  = workflow.workDir
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastp --version > v_fastp.txt
    """
}


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

/*
................................................................
                        PROCESSES
................................................................
*/

// Long-reads Quality Control

log.info" Quality control on Reads"

process longqc {
    /*
       longqc process to remove adapters and low quality sequences
    */
    publishDir "$params.outdir/quality_control/",
        mode: 'copy'
    tag "filter $sample_id"
    //memory '24 GB'
    
    input:
    tuple  val(sample_id), file(datasetFile)
    //val  threads

    output:
    tuple  val(sample_id), file("longqc_out/${sample_id}_trimmed.fastq.gz"), emit: lreads
    path("longqc_out/${sample_id}/${sample_id}.longqc.json"), emit: json
    path("longqc_out/${sample_id}/web_summary_${sample_id}.html")

    script:
    """
    mkdir longqc_out
    longQC.py sampleqc -x ont-rapid --sample_name ${sample_id} -o longqc_out/${sample_id} ${datasetFile} --trim_output longqc_out/${sample_id}_trimmed.fastq.gz --ncpu ${task.cpus}

    sed 's/QC_vals_longQC_sampleqc_//g' longqc_out/${sample_id}/QC_vals_longQC_sampleqc_${sample_id}.json > ${sample_id}.longqc.json
    mv ${sample_id}.longqc.json longqc_out/${sample_id}/
    """
}


// Short-reads Quality Control

process fastp {
    /*
       fastp process to remove adapters and low quality sequences
    */
    tag "filter $sample_id"
    publishDir "$params.outdir/quality_control/fastp_out",
        mode: 'copy'
    input:
    tuple val(sample_id), path(sreads_ch)

    output:
    tuple val(sample_id), path("${sample_id}_filt_*.fastq.gz"), emit: sreads
    path("${sample_id}.fastp.json"), emit: json


    script:
    """
    fastp -i ${sreads_ch[0]} -I ${sreads_ch[1]} \
      -o ${sample_id}_filt_1.fastq.gz -O ${sample_id}_filt_2.fastq.gz \
      --detect_adapter_for_pe -w ${task.cpus} -j ${sample_id}.fastp.json

    """
}

// Hybrid Assemblies
// --gsize 3.2M
log.info" Hybrid Genomels Assembly using $params.asb_tool "

if (params.asb_tool == 'dragonflye') {
    if (params.genomeSize == 0){
        log.error "No genome size specified. Necessary for Dragonflye (like 3.7M) assembly workflow"
        exit 1
    }
    process hb_dragonflye {
        tag { sample_id }
        publishDir "$params.outdir/assemblies/dflye",
            mode: 'copy'
        //module 'dragonflye/1.0.12'
	//cpus 8
       //memory '16 GB'

        input:
	//tuple val(sample_id), path(lreads)
	//tuple val(sample_id), path(sreads)
        tuple val(sample_id), path(join_l), path(lr)

        output:
        tuple val(sample_id), path("dflye_out/${sample_id}.fasta")

        //when:
	//isMode(['dragonflye'])

        script:
	"""
	mkdir -p dflye_out
        dragonflye --reads ${lr} --gsize ${params.genomeSize} --outdir ${sample_id} --polypolish 2 --R1 ${join_l[0]} --R2 ${join_l[1]} --cpus ${task.cpus} --force
	mv ${sample_id}/contigs.fa dflye_out/${sample_id}.fasta
        """
    }

}

// Create assembly with Canu

if (params.asb_tool == 'canu') {
    if (params.genomeSize == 0){
        log.error "No genome size specified. Necessary for Canu assembly workflow"
        exit 1
    }
    
    process hb_canu {
        tag {sample_id}
        publishDir "$params.outdir/assemblies/canu_asbl_out", mode: 'copy'
        //module 'canu/2.1:gnuplot/5.4.1'
        //module 'apptainer/latest'     
        input:
        tuple  val(sample_id), file(datasetFile)

        output:
        tuple val(sample_id), path("assembly_result_canu/${sample_id}.contigs.fasta")
        tuple val(sample_id), path("canu_results/*")
        
        script:
        """
        mkdir -p assembly_result_canu canu_results
        canu -p ${sample_id} -d ${sample_id} genomeSize=$params.genomeSize -nanopore ${datasetFile} 
        mv ${sample_id}/${sample_id}.contigs.fasta assembly_result_canu/
        mv ${sample_id} canu_results/
        """
        }
    
        
    // Map short reads to canu assembly with minimap2
    
    log.info" Alignment for Assembly polishing with Pilon "
    process minimap {
        tag { sample_id }
        publishDir "${params.outdir}/assemblies/minimap_out", mode: 'copy'
        
        input:
        tuple val(sample_id), path(asb_can), path(min_join) 
        //tuple val(sample_id), path(min_join)

        output:
        tuple val(sample_id), path("minimap_alignment_results/*")
        //tuple val(sample_id), path("short_reads_mapped_bam/*.sorted.bam")

        script:
        """
        mkdir -p minimap_alignment_results
        minimap2 -ax sr ${asb_can} ${min_join[0]} ${min_join[1]} > minimap_alignment_results/${sample_id}_aln.sam
        """
    }
    // bwa alignment - index and genome in the same direectory
    // ${sample_id}_index =`find -L ./ -name "${sample_id}*.amb" | sed 's/.amb//g'`
    // mv ${sample_id}.contigs* ${sample_id}_indx/
    process bwa_aln {
        tag { sample_id }
        publishDir "${params.outdir}/assemblies/bwa_out", mode: 'copy'
        input:
        tuple val(sample_id), path(asb_can), path(min_join)
        output:
        //tuple val(sample_id),path("*") 
        tuple val(sample_id), path("bwa_sam/${sample_id}_aln-pe.sam")
        script:
        """
        mkdir -p bwa_sam 
        bwa index $asb_can 
        bwa mem -t ${task.cpus} $asb_can ${min_join[0]} ${min_join[1]} > bwa_sam/${sample_id}_aln-pe.sam  
        """

    }
    // generating bam files, sorting and indexing
    
    process samtools_conv {

        tag { sample_id }
        publishDir "${params.outdir}/assemblies/samt_out", mode: 'copy'
        input:
        tuple val(sample_id), path(sam_aln)
        output:
        tuple val(sample_id), path("sr_mappedSorted_bam/*.sorted.bam")
        tuple val(sample_id), path("sr_indexed_bam/*.bam.bai")       
        
        script:
        """
        mkdir -p sr_mappedSorted_bam sr_indexed_bam 
        samtools view -h -b ${sam_aln} > ${sample_id}_aln.bam
        samtools sort ${sample_id}_aln.bam > sr_mappedSorted_bam/${sample_id}_aln.sorted.bam
        samtools index sr_mappedSorted_bam/${sample_id}_aln.sorted.bam    
        mv sr_mappedSorted_bam/*.bam.bai sr_indexed_bam/
        """
    } 

    // Polish assembly with pilon using assembly fasta and alignment
    // cp {${bam_indx},${alnbam}} bam_sortInd/
    process pilon {
        tag { sample_id }
        publishDir "${params.outdir}/assemblies/canu_asbl_out/pilon_out", mode: 'copy'

        input:
        tuple val(sample_id), path(can_asb), path(bam_indx), path(alnbam)
        //tuple val(sample_id), path(asb) //path(alnbam)
        //tuple val(sample_id),path(bam_indx)
        //tuple val(sample_id),path(alnbam)
        output:
        tuple val(sample_id), path("pilon_asb_scaffolds/*")
        // java -jar -Xmx96G pilon-1.24.jar
        script:
        """
        mkdir -p pilon_asb_scaffolds
        java -jar -Xmx16G /pilon/pilon.jar --genome ${can_asb} --bam ${alnbam} --outdir ${sample_id} 
        mv ${sample_id}/pilon.fasta ${sample_id}/${sample_id}.fasta
        mv ${sample_id}/${sample_id}.fasta pilon_asb_scaffolds/
        """

    }
    
}

// Collecting either assembly for downstream analyses
    log.info" Collecting Assemblies from $params.asb_tool for Downstream Analyses "
// rm -r /path/to/dir/*

process hbAssembly {
    publishDir "$params.outdir/assemblies/", mode: 'copy'
    tag { sample_id }
    input:
    tuple val(sample_id), path(whateverAsb)
    output:
    tuple val(sample_id), path("either_assembly/*")
    script:
    """
    mkdir -p either_assembly
    rm -r either_assembly;mkdir -p either_assembly
    cp -R ${whateverAsb} either_assembly/
    ls either_assembly
    """

    }


// Hybrid Unicycler

//if (params.asb_tool == 'unicycler') {

process hb_unicycler {

    tag { sample_id }
    
    publishDir "$params.outdir/assemblies/unicycler_out",
        mode: 'copy'

    input:
    tuple val(sample_id), path(join_l), path(lr)
    //tuple val(sample_id), path(lreads)
    //path(hb)

    output:
    tuple val("${sample_id}"), path("${sample_id}_assembly") // curly and $ added on 11/18/22
    tuple val("${sample_id}"), path("unic_asb/${sample_id}.fasta")
    script:
    """
    mkdir -p unic_asb
    unicycler -1 ${join_l[0]} -2 ${join_l[1]} -l ${lr} -o ${sample_id}_assembly --min_fasta_length 300 --keep 1 -t ${task.cpus} --min_kmer_frac 0.3 --max_kmer_frac 0.9 --verbosity 2
    mv ${sample_id}_assembly/assembly.fasta unic_asb/${sample_id}.fasta
    """
}

// Quast for any assembly
    log.info" Assembly Quality Metrics"
process quast {
    publishDir "${params.outdir}/quast_out", mode: 'copy'

    input:
    path(scaffolds)
    
    output:
    path("quast_results_out/*") // just one table

    script:
    """
    quast.py $scaffolds -o quast_results_out
    """
 }

// Species Identification

log.info" Species Identification"
process krakenMax {
    publishDir "$params.outdir/species-identification/kraken2Max_out",
    mode: 'copy'

    // show in the log which input file is analyzed
    tag { sample_id }
    input:
   //tuple val(sample_id), path(hbassembly)
    tuple val(sample_id), path(any_asb)

    output:
    tuple val(sample_id), path("kraken_report/${sample_id}.report")
    tuple val(sample_id), path("kraken_out/${sample_id}_kraken.out")

    script:
        if(!params.kradb && !params.krdbName){
			 """
        mkdir -p kraken_report kraken_out
        kraken2 --db /kraken2-db/minikraken2_v1_8GB/ --threads ${task.cpus} --use-names --report ${sample_id}.report --output ${sample_id}_kraken.out ${any_asb}
        mv ${sample_id}.report kraken_report/
        mv ${sample_id}_kraken.out kraken_out/
        """
	}
	else {
	"""
        mkdir -p kraken_report kraken_out
        kraken2 --db $params.kradb/$params.krdbName --threads ${task.cpus} --use-names --report ${sample_id}.report --output ${sample_id}_kraken.out ${any_asb}
        mv ${sample_id}.report kraken_report/
        mv ${sample_id}_kraken.out kraken_out/
        """
	}
	
}


process mash {
    publishDir "$params.outdir/species-identification/mash_out",
    mode: 'copy'
    tag { sample_id }
    input:
    tuple val(sample_id), path(any_asb)
    output:
    tuple val(sample_id), path("${sample_id}_distances_top10.tab")

    script:
    """
    mash sketch ${any_asb} -o ${sample_id}_sketch
    mash dist /db/RefSeqSketchesDefaults.msh ${sample_id}_sketch.msh > ${sample_id}_distances.tab
    sort -gk3 ${sample_id}_distances.tab | head > ${sample_id}_distances_top10.tab
    """
}

// Summary Reports
log.info" Kaiju, Mash, and Kraken Summary Reports"
// One direction for kaiju all Kaiju report files

// Adding all Kaiju report into a single directory
// Mash Summary report

// Adding all all report into a single directory
// Mash Summary report

process mash_summary{

publishDir "$params.outdir/species-identification/mashSummary_out",
        mode: 'copy'
    input:
    path(mash_grouped)

    output:
    path("*.txt")
    path("*.csv")

    script:
    """
    mash_summary_v02.py --mreppath ${mash_grouped}
    """
}

// Adding all kraken reports into one directory

// Kraken Summmary directory of reports as input

process kraken_summary{

    publishDir "$params.outdir/species-identification/krakenSummary_out",
        mode: 'copy'
    input:
    path(make_dir_4kreport)

    output:
    path("krakenMax_report.txt")
    path("sample_gen_tuple.txt"), emit: sample_gen
    path("tp_sampTaxid.txt"), emit:samp_taxId

    script:
    """
    kraken_summaryv2_5tuples.py --reppath ${make_dir_4kreport} --out krakenMax_report.txt --tupleGfile sample_gen_tuple.txt --taxonId tp_sampTaxid.txt
    """
}

// Genome Annotation

log.info" Genome Annotation"

process prokka {
    tag { sample_id }
    
    publishDir "$params.outdir/annotation/prokka_out", 
        mode: 'copy'
    cpus 8
    memory '16 GB'
    module 'prokka' 
    input:
    tuple val(sample_id), path(any_asb)      
    
    output:
    //tuple val(sample_id), path("${sample_id}/report.{pdf,html}")
    path("${sample_id}")

    script:
    """
  
    prokka --cpus ${task.cpus} --fast --outdir ${sample_id} --prefix ${sample_id} --force --compliant ${any_asb} --
    """
}
// Adding all prokka results into one dir first
process prokka_summary {
    publishDir "$params.outdir/annotation/prokka_summary", mode: 'copy'

    input:
        path(tx)

    output:
        path("*")

    script:
    """
    prokka_summarygbk_v01.py --pkanres ${tx}
    """
}


log.info" Genome Annotation Summary Reports from Prokka, Bakta, and Pgap"

// Creating ymal metadata files by sample using sampleID_species tuples from kraken_summary

process metaDt_gener {
    
    publishDir "$params.outdir/annotation/metaDt_dir",
        mode: 'copy'
    input:    
        path(sample_gen)
        path(meta)
    
    output:
        path("*")
    
    script:
    """
    ymal_editbySample_v02.1.py --sample_species ${sample_gen} --input_yaml ${meta[1]} --input_submol ${meta[0]} --input_ext fasta    
    """
}

// grouping all the metadata directories into one directory (parent->child)
process grouping_metaDtDir{
    publishDir "$params.outdir/annotation/",
        mode: 'copy'
    debug true
    input:
        path(yp)
    output:
        path("dir_metaDtDir")
    script:
    """
    echo "[make_dir]"
    mkdir -p dir_metaDtDir
    
    mv ${yp} dir_metaDtDir/
    echo "[user_dir]"
    ls dir_metaDtDir
    """
}
// Having all asemblies into a single dir for copying into metaDt_dir
// Copying assembly files from assbl_dir into same directories as yaml metada

process groupingfor_pgap {
    
    //tag {sample_id}
    publishDir "$params.outdir/annotation",
        mode: 'copy'
    input:
        //tuple val(sample_id),path(asb)
        //path(asb)
        path(assbl_dir)
        path(dir_metdir)
        //path(mtd)

    output:
        //tuple val(sample_id),path("${sample_id}")
        path("$dir_metdir"), emit: asb_met_dir // distination directory
        
    script:
    """
    grouping_asb_met04.py --asb_dir $assbl_dir --met_dir $dir_metdir
    ls $dir_metdir        
    """
}

// NCBI pgap annotation // May be using the direct command

process pgap_antn {
    //tag {sample_id}
    publishDir "$params.outdir/annotation/pgap_out",
        mode: 'copy'
    input:
        //tuple val(sample_id), path(wp)
        path(asb_met_dir)
    output:
        //tuple val(sample_id), path("${sample_id}_pgap/*")
        path("*")
    script:
    """
    pgap_run_v3.py --ppgroup_dir $asb_met_dir
    """
}

// pgap Summary Report

// Grouping pgap result into one directory first

process pgap_summary{
publishDir "$params.outdir/annotation/pgapSummary_out",
        mode: 'copy'
    input:
    path(pg)

    output:
    path("*.txt")
    path("*.csv")
    path("*.xlsx")

    script:
    """
    pgap_summary_v01.py --pgapres ${pg}
    """
}

// Average Nucleotide Identity
log.info" Average Nucleotide Identity "

process ani {
    //tag {sample_id}
    publishDir "$params.outdir/ani_out/ani_final", mode: 'copy'
    
    input:
        path(asb_met_dir)
        
    output:
        path("*")
        
    script:
    """
    pgap_tax_v2.py --ppgroup_dir ${asb_met_dir}
    """
}

// Adding all ani report to a single directory first
process ani_summary {
    publishDir "$params.outdir/ani_out/ani_final", mode: 'copy'

    input:
        path(tx)

    output:
        path("*")

    script:
    """
    ani_summary_v02.py --aniReppath ${tx}
    """
}

// Antimicromial resistance

log.info" Detection of Antimicrobial Resistance Genes"

process amrfinderplus {
    tag { sample_id }
    publishDir "$params.outdir/amrfinder/", mode: 'copy'
    input:
    tuple val(sample_id), path(any_asb)
    output:
    tuple val(sample_id), path("${sample_id}.tsv")

     script:
    """
    amrfinder --nucleotide ${any_asb} --plus --threads ${task.cpus} -o ${sample_id}.tsv
    """
 }

// Importing Modules
//include {grouping_intodir as grouping_intodir_kaiju } from './flisochar_module_v02.nf'
include {grouping_intodir as grouping_intodir_mash} from './flisochar_module_v02.nf'
include {grouping_intodir as grouping_intodir_kraken} from './flisochar_module_v02.nf'
//include {grouping_intodir as grouping_intodir_metaDt} from './flisochar_module_v02.nf'
include {grouping_intodir as grouping_intodir_prokka} from './flisochar_module_v02.nf'
//include {grouping_intodir as grouping_intodir_bakta} from './flisochar_module_v02.nf'
include {grouping_intodir as grouping_intodir_asbl} from './flisochar_module_v02.nf'
include {grouping_intodir as grouping_intodir_pgap} from './flisochar_module_v02.nf'
include {grouping_intodir as grouping_intodir_ani} from './flisochar_module_v02.nf'

// Workflow

workflow ontflow {
    // Quality control
    longqc( lreads_ch)
    fastp( sreads_ch )
    lqc_ch = longqc.out.lreads
    sqc_ch = fastp.out.sreads
    lsch_join = sqc_ch.join(lqc_ch) // long reads and short reads channel
    // Dragonflye Assembly
    if (params.asb_tool == 'dragonflye') {
        df_ch = hb_dragonflye(lsch_join)
        hbCollect_ch = hbAssembly(df_ch)
        quast_ch = quast(hbCollect_ch.map{it -> it[1]}.collect())
        
        krakenM_ch = krakenMax(hbCollect_ch) // Kraken species ID
        
        //kaiju_ch = kaiju(hbCollect_ch)
        //kaijuRep_list = kaiju_ch.map{it -> it[1]}.collect()
        //kaiju_grouped_ch = grouping_intodir_kaiju(kaijuRep_list)
        //kaiju_summary_ch = kaiju_summary(kaiju_grouped_ch)
      
        mash_ch = mash(hbCollect_ch)
        mashRes_list = mash_ch.map{it -> it[1]}.collect()
        mash_grouped_ch = grouping_intodir_mash(mashRes_list)
        mash_summary_ch = mash_summary(mash_grouped_ch)
        
        prokkaAnt_ch = prokka(hbCollect_ch)
        prokkaRes_list = prokkaAnt_ch.collect()
        prokka_grouped_ch = grouping_intodir_prokka(prokkaRes_list)
        prokka_summary_ch = prokka_summary(prokka_grouped_ch)

        //bakta_ch = bakta(hbCollect_ch)
        //baktaRes_list = bakta_ch.collect()
        //bakta_grouped_ch = grouping_intodir_bakta(baktaRes_list)
        //bakta_summary_ch = bakta_summary(bakta_grouped_ch)

      	amrfinderplus_ch = amrfinderplus(hbCollect_ch)
        
        krakenRep_list = krakenM_ch[0].map{it -> it[1]}.collect()
        kraken_grouped_ch = grouping_intodir_kraken(krakenRep_list)
        kraken_summary_ch = kraken_summary(kraken_grouped_ch)
        
        meta_gener_ch = metaDt_gener(kraken_summary_ch[1], meta_ch)
        metDir_list = meta_gener_ch.collect()
        metdir_grouped_ch = grouping_metaDtDir(metDir_list)
        
        asbly_list = hbCollect_ch.map{it -> it[1]}.collect()
        asbly_grouped_ch = grouping_intodir_asbl(asbly_list)
        
        groupingfor_pgap_ch = groupingfor_pgap(asbly_grouped_ch, metdir_grouped_ch)
        pgap_antn_ch = pgap_antn(groupingfor_pgap_ch)
        pgap_resList = pgap_antn_ch.collect()
        pgap_resgrouped_ch = grouping_intodir_pgap(pgap_resList)
        pgap_summary_ch = pgap_summary(pgap_resgrouped_ch)
        
        ani_ch = ani(groupingfor_pgap_ch)
        aniRes_list = ani_ch.collect()
        ani_grouped_ch = grouping_intodir_ani(aniRes_list)
        ani_summary_ch = ani_summary(ani_grouped_ch)
    }
    
    // Unicycler Assembly
    if (params.asb_tool == 'unicycler') {
        lsch_join = sqc_ch.join(lqc_ch) // long reads and short reads channel
        unic_ch = hb_unicycler(lsch_join)
        hbCollect_ch = hbAssembly(unic_ch[1])
        quast_ch = quast(hbCollect_ch.map{it -> it[1]}.collect()) // asb metrics
        
        krakenM_ch = krakenMax(hbCollect_ch) // Kraken species ID
                       
        mash_ch = mash(hbCollect_ch)
        mashRes_list = mash_ch.map{it -> it[1]}.collect()
        mash_grouped_ch = grouping_intodir_mash(mashRes_list)
        mash_summary_ch = mash_summary(mash_grouped_ch)
        
        prokkaAnt_ch = prokka(hbCollect_ch)
        prokkaRes_list = prokkaAnt_ch.collect()
        prokka_grouped_ch = grouping_intodir_prokka(prokkaRes_list)
        prokka_summary_ch = prokka_summary(prokka_grouped_ch) 
        
        amrfinderplus_ch = amrfinderplus(hbCollect_ch)
        
        krakenRep_list = krakenM_ch[0].map{it -> it[1]}.collect()
        // make_dir_4kreport_ch = make_dir_4kreport(list_kreport)
        kraken_grouped_ch = grouping_intodir_kraken(krakenRep_list)
        kraken_summary_ch = kraken_summary(kraken_grouped_ch)
        
        meta_gener_ch = metaDt_gener(kraken_summary_ch[1], meta_ch)
        metDir_list = meta_gener_ch.collect()
        metdir_grouped_ch = grouping_metaDtDir(metDir_list)
        
        asbly_list = hbCollect_ch.map{it -> it[1]}.collect()
        asbly_grouped_ch = grouping_intodir_asbl(asbly_list)
        groupingfor_pgap_ch = groupingfor_pgap(asbly_grouped_ch, metdir_grouped_ch)
        
        pgap_antn_ch = pgap_antn(groupingfor_pgap_ch)
        pgap_resList = pgap_antn_ch.collect()
        pgap_resgrouped_ch = grouping_intodir_pgap(pgap_resList)
        pgap_summary_ch = pgap_summary(pgap_resgrouped_ch)
        
        ani_ch = ani(groupingfor_pgap_ch)
        aniRes_list = ani_ch.collect()
        ani_grouped_ch = grouping_intodir_ani(aniRes_list)
        ani_summary_ch = ani_summary(ani_grouped_ch)
    }   

    // Canu Assembly
    if (params.asb_tool == 'canu') {
        ca_ch = hb_canu(lreads_ch)
        //hbAssembly(ca_ch[0])
        
        minim_join_ch = ca_ch[0].join(sreads_ch) // joining canu asb ch and sr ch for alignment
        minim_ch = minimap(minim_join_ch) // minmap2 alignment
        //samtMi_ch = samtools_conv(minim_ch)
        bwa_ch = bwa_aln(minim_join_ch) // bwa alignment
        samtBw_ch =  samtools_conv(bwa_ch) // sam processing aln
        
        // Polishing from BWA alignment 
        pilonBw_join1 = ca_ch[0].join(samtBw_ch[1]) // join asb and idx
        pilonBw_join2 = pilonBw_join1.join(samtBw_ch[0]) // join asb,indx-bam, and sorted_bam
        pilonBw_ch = pilon(pilonBw_join2) // Polishing 
        
        // Collecting the polished assembly
        hbCollect_ch = hbAssembly(pilonBw_ch)
        
        //Polishing from  Minimap2 alignment 
        //pil_join_ch = ca_ch[0].join(samtMi_ch[1]) // join lr asbl and minimap2 alignment
        //pil_join_ch.view()
        //pil2_join_ch=pil_join_ch.join(samt_ch[0])
        //pilon_polish_ch = pilon(pil2_join_ch)

        // Assembly quality assessment
        quast_ch = quast(hbCollect_ch.map{it -> it[1]}.collect())
        krakenM_ch = krakenMax(hbCollect_ch) // Kraken species ID
        
        prokkaAnt_ch = prokka(hbCollect_ch)
        prokkaRes_list = prokkaAnt_ch.collect()
        prokka_grouped_ch = grouping_intodir_prokka(prokkaRes_list)
        prokka_summary_ch = prokka_summary(prokka_grouped_ch)

        mash_ch = mash(hbCollect_ch)
        mashRes_list = mash_ch.map{it -> it[1]}.collect()
        mash_grouped_ch = grouping_intodir_mash(mashRes_list)
        mash_summary_ch = mash_summary(mash_grouped_ch)
        
        
        amrfinderplus_ch = amrfinderplus(hbCollect_ch)
        
        krakenRep_list = krakenM_ch[0].map{it -> it[1]}.collect()
        kraken_grouped_ch = grouping_intodir_kraken(krakenRep_list)
        kraken_summary_ch = kraken_summary(kraken_grouped_ch)
        
        meta_gener_ch = metaDt_gener(kraken_summary_ch[1], meta_ch)
        metDir_list = meta_gener_ch.collect()
        metdir_grouped_ch = grouping_metaDtDir(metDir_list)
        
        asbly_list = hbCollect_ch.map{it -> it[1]}.collect()
        asbly_grouped_ch = grouping_intodir_asbl(asbly_list)    
        groupingfor_pgap_ch = groupingfor_pgap(asbly_grouped_ch, metdir_grouped_ch)
        
        pgap_antn_ch = pgap_antn(groupingfor_pgap_ch)
        pgap_resList = pgap_antn_ch.collect()
        pgap_resgrouped_ch = grouping_intodir_pgap(pgap_resList)
        pgap_summary_ch = pgap_summary(pgap_resgrouped_ch)
        
        ani_ch = ani(groupingfor_pgap_ch)
        aniRes_list = ani_ch.collect()
        ani_grouped_ch = grouping_intodir_ani(aniRes_list)
        ani_summary_ch = ani_summary(ani_grouped_ch)

    }
}   

/*
 * Completion email notification
 */

workflow.onComplete {
    log.info "[flisochar] Pipeline Complete! visit the following directory in your current directory --> $params.outdir/"

}
// email notification

/*
   def msg = """\
            Pipeline execution summary
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            """
            .stripIndent()
    sendMail(to: $params.email, subject: 'Flisochar execution', body: msg)
*/