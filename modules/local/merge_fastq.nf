/*
 * Define "merge_fastq" process
 *
 */
process MERGE_FASTQ {
    tag "merge_fastq"
    publishDir params.sampledir,  mode: "copy", pattern: "merged_samples.csv"
    publishDir params.logdir, mode: "copy", pattern: "merge_fastq.out"


    input:
    path (ch_samples)

    output:
    path("*.fastq.gz"),                 emit: reads
    path("merge_fastq.out"),            emit: log
    path("merged_samples.csv"),          emit: ch_merged_samples


    script:
    """
    sed -e 's/[[:space:]]*\$//' ${ch_samples} | sed 's/"*\$//g' | sed 's/^"*//g' | grep -v "^#" > samples_tmp.txt

    echo "sample,fastq_1,fastq_2,strandedness,groupname" > merged_samples.csv 

    if [ ! -d  ${params.merged_fastq} ]; then
        mkdir -p  ${params.merged_fastq}
    fi

    while read line; do
        if [[ ! "\$line" =~ "^#" ]]; then
            read group control replicate spike email r1s r2s <<< \$line    
            sample=\${group}_\${replicate}
            merged_r1=\${sample}_r1.fastq.gz
            merged_r2=\${sample}_r2.fastq.gz
            cat \${r1s//,/ } > ${params.merged_fastq}/\${merged_r1}
            cat \${r2s//,/ } > ${params.merged_fastq}/\${merged_r2}
            ln -s ${params.merged_fastq}/\${merged_r1} \${merged_r1}
            ln -s ${params.merged_fastq}/\${merged_r2} \${merged_r2}


            echo "\${sample},${params.merged_fastq}/\${merged_r1},${params.merged_fastq}/\${merged_r2},unstranded,\${group}"

        fi
    done < samples_tmp.txt  >> merged_samples.csv
   
 
    cat .command.log > merge_fastq.out
    """
}

