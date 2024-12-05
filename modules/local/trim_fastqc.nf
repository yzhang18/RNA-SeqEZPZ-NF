/*
 * Define "trim" process
 * Note trimming and fastqc are done by trim_galore
 *
 */
process TRIM_FASTQC {
    tag "$meta.id"
    label TRIM_FASTQC
    publishDir params.trim,  mode: "copy", pattern: "*.{fq.gz,txt}"
    publishDir params.fastqc, mode: "copy", pattern: "*.{html,zip}"
    publishDir params.logdir, mode: "copy", pattern: "trim_fastqc_*.out"

    input:
    tuple val(meta), path(reads)
    path (chr)

    output:
    tuple val(meta), path("*.fq.gz"),        emit: reads
    //path("*.fq.gz"),        emit: reads
    tuple val(meta), path("*report.txt"),    emit: report
    tuple val(meta), path("*.zip"),          emit: zip
    tuple val(meta), path("*.html"),         emit: html
    path("trim_fastqc_*.out"),            emit: log


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -x
    source activate rna_env
    trim_galore --version

    #ls -ld ${params.outdir}
    trim_galore -q 20 --length 20 --illumina --cores ${params.ncpus_trim} --paired ${reads[0]} ${reads[1]} --fastqc

    conda deactivate
    set +x

    cat .command.log > trim_fastqc_${prefix}.out
    """
}

