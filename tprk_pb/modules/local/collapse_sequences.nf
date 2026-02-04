process COLLAPSE_SEQUENCES {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/fedora/python-312:312'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("${meta.id}.filtQ*.fasta"), emit: fastas
    tuple val(meta), path ("${meta.id}.collapse.tsv"), emit: summary

    script:
    """
    dereplicate.py ${meta.id}.fw.trim.filtQ20.fastq ${meta.id}.filtQ20.fasta 20 > ${meta.id}.filtQ20.collapse.tsv
    dereplicate.py ${meta.id}.fw.trim.filtQ30.fastq ${meta.id}.filtQ30.fasta 30 > ${meta.id}.filtQ30.collapse.tsv
    dereplicate.py ${meta.id}.fw.trim.filtQ40.fastq ${meta.id}.filtQ40.fasta 40 > ${meta.id}.filtQ40.collapse.tsv
    
    paste -d'\t' <(head -n1 ${meta.id}.filtQ20.collapse.tsv) <(head -n1 ${meta.id}.filtQ30.collapse.tsv) <(head -n1 ${meta.id}.filtQ40.collapse.tsv) > ${meta.id}.collapse.tsv
    paste -d'\t' <(tail -n1 ${meta.id}.filtQ20.collapse.tsv) <(tail -n1 ${meta.id}.filtQ30.collapse.tsv) <(tail -n1 ${meta.id}.filtQ40.collapse.tsv) >> ${meta.id}.collapse.tsv

    """
}