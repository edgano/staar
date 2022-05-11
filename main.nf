#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/
========================================================================================

----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

//params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

//WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

// include { ASSOCIATIONS } from "${projectDir}/workflows/associations"
// include { STAAR } from "${projectDir}/workflows/staar"
include { TEST } from "${projectDir}/workflows/test"
//
// WORKFLOW: Run main nf-core/associations analysis pipeline
//
workflow NFCORE_ASSOCIATIONS {
    //ASSOCIATIONS ()
    TEST()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
/*
workflow {
    NFCORE_ASSOCIATIONS ()
}
*/
/*
========================================================================================
    THE END
========================================================================================
*/
