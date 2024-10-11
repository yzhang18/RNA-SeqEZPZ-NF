/*
 * Define "multiqc" process
 *
 *
 */
process MULTIQC1 {
    publishDir params.logdir, mode: "copy", pattern: "multiqc1.log"

    input:
    path("*.log")

    output:
    path ("multiqc1.log"),        emit: log
    path ("readlength.txt"),      emit: readlength

    script:
    """
    source activate rna_env
    export LANG=en_US.UTF-8 && export LC_ALL=en_US.UTF-8
    log_dir=\$(pwd)
    #cd ${params.multiqc}
    cd /mnt/outputs/fastqc_rslt
    multiqc --version
    multiqc -f .

    readlength=\$(awk 'BEGIN{n=0} {sum += \$15; n++} END {if (n>0) print int(sum/(n-1));}' \
         ./multiqc_data/multiqc_fastqc.txt)
    echo \$readlength > \${log_dir}/readlength.txt
    conda deactivate

    cd \${log_dir}
    cat .command.log > multiqc1.log
    """
}

