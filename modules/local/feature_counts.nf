/*
 *
 * define process feature_counts
 *
 */
process FEATURE_COUNTS {
    tag "$meta.id"
    label "hi_cpus"
    publishDir params.star_pass2, mode: "copy", pattern: "*{counts.txt,counts.txt.summary}"
    publishDir params.logdir, mode: "copy", pattern: "feature_counts_*.log"

    input:
    tuple val(meta), path(bam_sorted_file)
    val (gtf)

    output:
    tuple val(meta), path("*_counts.txt")        , emit: counts
    tuple val(meta), path("*_counts.txt.summary"), emit: counts_summary
    path("feature_counts_*.log")                 , emit: log

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    """
    source activate rna_env
    featureCounts -v
    set -x
    featureCounts -p \
            -a /ref/${params.genome}/${gtf} \
            -t exon \
            -g gene_id \
            -T ${task.cpus} \
            -o ${prefix}_counts.txt ${bam_sorted_file}

    set +x
    conda deactivate

    cat .command.log > feature_counts_${prefix}.log
    """
}

