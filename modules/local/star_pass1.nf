/*
 * Define "star_pass1" process
 *
 *
 */
process STAR_PASS1 {
    tag "$meta.id"
    label "hi_mem_cpus"
    publishDir params.logdir, mode: "copy", pattern: "star_pass1_*.log"


    input:
    tuple val(meta), path(reads)
    path (index)


    output:
    tuple val(meta), path("*d.out.sam")       , emit: sam
    tuple val(meta), path("*Log.final.out")   , emit: log_final
    tuple val(meta), path("*Log.out")         , emit: log_out
    tuple val(meta), path("*Log.progress.out"), emit: log_progress
    tuple val(meta), path("*.tab")            , optional:true, emit: tab
    path("*.sjdb")                            , optional:true, emit: sjdb
    path("star_pass1_*.log")                  , emit: log

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    source activate rna_env
    set -x 
    STAR --genomeDir ${index} \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${prefix}_
        
    genome_for_pass2.sh ${prefix}
    set +x
    conda deactivate

    cat .command.log > star_pass1_${prefix}.log
    """
}   


