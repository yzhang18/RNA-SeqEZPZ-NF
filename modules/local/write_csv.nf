/*
 *
 * define process WRITE_CSV
 *
 */
process WRITE_CSV {
    publishDir "$params.sampledir", mode: "copy", pattern: "*.csv"

    input:
    val (items)

    output:
    path ("samples.csv"), emit: samples_csv

    // Write sample information to a csv file
    script:
    """ 
    # remove the first and last characters "[" and "]"
    sample_info=\$(echo -n "${items}" | tail -c +2 | head -c -1)
 
    echo "sample,fastq_1,fastq_2,strandedness,groupname,controlname,replicatename" > samples.csv
    export IFS=","
    set -- \${sample_info}

    while [ \$# -ne 0 ]  
    do
        echo "\$1,\$2,\$3,\$4,\$5,\$6,\$7" >> samples.csv
        shift 7
    done
    """ 
}

