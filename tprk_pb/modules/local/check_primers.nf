process CHECK_PRIMERS {
    tag "$meta.id"
    label 'process_single'
    container 'biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(meta), path(raw_fastq), path(fw_fastq), path(trim_fastq)


    output:
    tuple val(meta), path ("${meta.id}.primers_summary.tsv"), emit: summary

    script:
    """
    check_primers.sh ${raw_fastq} ${fw_fastq} ${trim_fastq} ${params.primer_error} > ${meta.id}.primers_summary.tsv

    """
}