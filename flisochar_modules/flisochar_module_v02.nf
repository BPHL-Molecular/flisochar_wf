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

// Input parameters 

nextflow.enable.dsl = 2

//params.rep = ""
//params.outdir="$baseDir/flisoch_out"

// Channels

//rep_ch = Channel
//        .fromPath(params.rep, checkIfExists: true)

// Processes
process grouping_intodir{
    //publishDir "$params.outdir/",
        //mode: 'copy'
    debug true
    input:
    path(x)
    output:
    path("grouping_dir")
    script:
    """
    echo "[make_dir]"
    mkdir -p grouping_dir
    
    mv ${x} grouping_dir/
    echo "[user_dir]"
    ls grouping_dir
    """
}

//workflow{
//list_rep_ch = rep_ch.collect() // in case of dir ch
//list_rep_ch = rep_ch.map{it -> it[1]}.collect() // in case of files ch
//grouping_intodir_ch = grouping_intodir(list_rep_ch)
//rep_ch.collect().view()
//}


//workflow test1{
//list_rep_ch = rep_ch.collect()
//grouping_intodir_ch1 = grouping_intodir(list_rep_ch)

//}