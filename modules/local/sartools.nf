/*
 *
 * define process sartools
 *
 */
process SARTOOLS {
    tag ""
    label "very_hi_mem"
    publishDir params.logdir, mode: "copy", pattern: "sartools.log"

    input:
    path("*_counts.txt")

    output:
    path("sartools.log"), emit: log 


    script:
    """
    set +eu
    set -x
    source activate sartools_env
    log_dir=\$(pwd)
    cd /mnt/outputs

    Rscript /scripts/template_script_DESeq2_simg.R padj=${params.padj} email=${params.email} batch_adjust=${params.batch_adjust} -e "q(status=1)"
    if [[ \$? -eq 1 ]];then
        exit 1
    fi

    conda deactivate
    set +x
    #set -eu

    cd \${log_dir}
    cat .command.log > sartools.log
    """
}


