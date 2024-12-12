/*
 * define star_pass2 process
 *
 */
process STAR_PASS2 {
    tag "$meta.id"
    label "hi_mem_cpus"
    publishDir params.star_pass2, mode: "copy", pattern: "*{out,tab,txt,summary,STARgenome}"
    publishDir params.bw_files, mode: "copy", pattern: "*.bw"
    publishDir params.logdir, mode: "copy", pattern: "star_pass2_*.out"

    input:
    tuple val(meta), path(reads)
    path (sjFiles)
    path (gsize)
    path (index)

    output:
    tuple val(meta), path("*Log.final.out")   , emit: log_final
    tuple val(meta), path("*Log.out")         , emit: log_out
    tuple val(meta), path("*Log.progress.out"), emit: log_progress
    tuple val(meta), path("*sortedByCoord.out.bam")  , optional:true, emit: bam_sorted
    tuple val(meta), path("*sortedByCoord.out.bam.bai"), emit: bai
    tuple val(meta), path("*.tab")                   , optional:true, emit: tab
    tuple val(meta), path("*__STARgenome")           , emit: stargenome
    tuple val(meta), path("*_norm.bw")               , emit: bw
    path("star_pass2_*.out")                         , emit: log



    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    source activate rna_env
    echo "STAR version " `STAR --version`
    set -x


    STAR --genomeDir ${index} \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate Unsorted \
        --sjdbFileChrStartEnd ${sjFiles} \
        --limitSjdbInsertNsj 4000000 \
        --outFileNamePrefix ${prefix}_
    set +x
    conda deactivate

    # create bam index and normalized bams
    source activate samtools_env
    samtools --version
    set -x
    samtools index ${prefix}_Aligned.sortedByCoord.out.bam
    set +x
    conda deactivate

    source activate deeptools_env
    bamCoverage --version
    set -x
    read genome_size < ${gsize}
    echo \$genome_size
    bamCoverage -b ${prefix}_Aligned.sortedByCoord.out.bam \
            -o ${prefix}_norm.bw \
            -p ${task.cpus} \
            --binSize 10 \
            --smoothLength 20 \
            --normalizeUsing BPM \
            --effectiveGenomeSize \${genome_size}
    set +x
    conda deactivate

    cat .command.log > star_pass2_${prefix}.out

    """
}

