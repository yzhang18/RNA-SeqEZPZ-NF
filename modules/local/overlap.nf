/*
 *
 * define process overlap
 *
 */
process OVERLAP {
    publishDir params.shiny_app,  mode: "copy", pattern: "*.{sh,R}"
    publishDir params.logdir, mode: "copy", pattern: "overlap.log"

    input:
    path("*.log")

    output:
    path("*.sh")
    path("*.R")
    path("overlap.log")

    script:
    """
    cp /scripts/app.R ./
    cp /scripts/app_simg.sh ./
    cp /scripts/run_overlap.sh ./

    cat .command.log > overlap.log
    """
}


