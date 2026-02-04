/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowTprkpb.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CHECK_PRIMERS } from '../modules/local/check_primers'
include { COLLAPSE_SEQUENCES } from '../modules/local/collapse_sequences'
include { FILTER_QUALITY } from '../modules/local/filter_quality'
include { ORIENT_VALIDATE } from '../modules/local/orient_validate'
include { SUMMARY } from '../modules/local/summary'
include { FINAL_SUMMARY } from '../modules/local/final_summary'
include { CUTADAPT } from '../modules/nf-core/cutadapt/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TPRKPB {

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )

    ORIENT_VALIDATE (
        INPUT_CHECK.out.reads
    )

    CUTADAPT (
        ORIENT_VALIDATE.out.fastq
    )

    INPUT_CHECK.out.reads
    .join(ORIENT_VALIDATE.out.fastq)
    .join(CUTADAPT.out.reads)
    .set{ ch_all_fastqs }

    CHECK_PRIMERS (
        ch_all_fastqs
    )

    FILTER_QUALITY (
        CUTADAPT.out.reads
    )

    COLLAPSE_SEQUENCES (
        FILTER_QUALITY.out.fastqs
    )

    CHECK_PRIMERS.out.summary
    .join(FILTER_QUALITY.out.summary)
    .join(COLLAPSE_SEQUENCES.out.summary)
    .set{ ch_all_summaries }

    SUMMARY (
        ch_all_summaries
    )

    FINAL_SUMMARY (
        SUMMARY.out.summary.collect()
    )
}


