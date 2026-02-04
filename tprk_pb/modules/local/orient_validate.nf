process ORIENT_VALIDATE {
    tag "$meta.id"
    label 'process_single'
    container 'biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(meta), path(reads)


    output:
    tuple val(meta), path ("${meta.id}.fw.fastq"), emit: fastq
    tuple val(meta), path ("${meta.id}.fw.rejected.fastq"), emit: rejected
    tuple val(meta), path ("${meta.id}.fw.log"), emit: log

    script:
    """
    flip_to_forward.sh ${reads[0]} ${meta.id}.fw.fastq ${params.primer_error} 1400 1700 > ${meta.id}.fw.log

    """
}