/*
 *
 * define combinebw process
 *
 */
process COMBINEBW {
    tag "$groupid"
    publishDir params.bw_files, mode:"copy", pattern: ""
    publishDir params.logdir, mode: "copy", pattern: "combinebw_*.log"

    input:
    path (chr)
    tuple val(groupid), path (bw_files)

    output:
    tuple val(groupid), path("*_comb.bw"), emit: comb_bw
    path("combinebw_*.log")             , emit: log

    script:
    """
    source activate wiggletools_env
    set -x
    wiggletools mean ${bw_files} | \
        wigToBigWig stdin ${chr} ${groupid}_comb.bw
    set +x
    conda deactivate

    cat .command.log > combinebw_${groupid}.log
    """
}

