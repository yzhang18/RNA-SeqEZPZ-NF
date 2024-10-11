/*
 * Define "multiqc" process
 *
 *
 */
process MULTIQC2 {
    publishDir params.logdir, mode: "copy", pattern: "multiqc2.log"

    input:
    path("*.log")

    output:
    path ("multiqc2.log"), emit: log

    script:
    """
    source activate rna_env
    export LANG=en_US.UTF-8 && export LC_ALL=en_US.UTF-8
    log_dir=\$(pwd)
    #cd ${params.multiqc}
    cd /mnt/outputs/fastqc_rslt
    multiqc --version
    multiqc -f ../
    conda deactivate

    cd \${log_dir}
    cat .command.log > multiqc2.log
    """
}

