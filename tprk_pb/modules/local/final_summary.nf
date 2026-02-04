process FINAL_SUMMARY {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'nf-core/ubuntu:20.04' }"

    input:
    path summary_tsv_files

    output:
    path "final_summary.tsv", emit: summary

    """
    awk 'FNR==1 && NR==1{print} FNR==2{print}' *.tsv > final_summary.tsv

    """ 
}