process SUMMARY {
    tag "$meta.id"
    label 'process_single'
    container 'biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(meta), path(primer_summary), path(filter_summary), path(collapse_summary)


    output:
    path "${meta.id}.summary.tsv", emit: summary

    script:
    """
    paste -d'\t' <(printf "sample") <(head -n1 ${primer_summary}) <(head -n1 ${filter_summary}) <(head -n1 ${collapse_summary}) > ${meta.id}.summary.tsv
    paste -d'\t' <(printf "${meta.id}") <(tail -n1 ${primer_summary}) <(tail -n1 ${filter_summary}) <(tail -n1 ${collapse_summary}) >> ${meta.id}.summary.tsv

    """    
}