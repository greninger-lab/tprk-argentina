process FILTER_QUALITY {
    tag "$meta.id"
    label 'process_single'
    container 'biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(meta), path(trim_fastq)


    output:
    tuple val(meta), path ("${meta.id}.fw.trim.filtQ*.fastq"), emit: fastqs
    tuple val(meta), path ("${meta.id}.filt.tsv"), emit: summary

    script:
    """
    filter_quality.sh ${trim_fastq} ${meta.id}.fw.trim.filtQ20.fastq 20 > ${meta.id}.filtQ20.tsv
    filter_quality.sh ${trim_fastq} ${meta.id}.fw.trim.filtQ30.fastq 30 > ${meta.id}.filtQ30.tsv
    filter_quality.sh ${trim_fastq} ${meta.id}.fw.trim.filtQ40.fastq 40 > ${meta.id}.filtQ40.tsv

    paste -d'\t' <(head -n1 ${meta.id}.filtQ20.tsv) <(head -n1 ${meta.id}.filtQ30.tsv) <(head -n1 ${meta.id}.filtQ40.tsv) > ${meta.id}.filt.tsv
    paste -d'\t' <(tail -n1 ${meta.id}.filtQ20.tsv) <(tail -n1 ${meta.id}.filtQ30.tsv) <(tail -n1 ${meta.id}.filtQ40.tsv) >> ${meta.id}.filt.tsv

    """
}