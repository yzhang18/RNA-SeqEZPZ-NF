/*
 *
 * define process get_sample
 *
 */
process GET_SAMPLE {
    label "short_time"
    publishDir "$params.outdir", mode: "copy", pattern: "*.txt"

    input:
    path (samples_txt)

    output:
    path ("tmp.csv"),  emit: samples
    path (samples_txt), emit: txt
    //path ("email.txt"), emit: email
    env (email), emit:email

    script:
    """
    read_samples.pl ${samples_txt} tmp.csv email.txt
    email=\$(cat email.txt)
    """
}

