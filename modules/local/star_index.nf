/*
 * Define the "index" process
 * using tools STAR 
 */
process STAR_INDEX {
    label "star"
    //publishDir params.star_index, mode: "copy", pattern: "STAR_index"

    input:
    val (fasta) 
    val (gtf) 

    output:
    path ("genome.chrom.sizes"), emit: chr
    path ("STAR_index"), emit: index

    script:
    """

    set -x

    # If STAR_index/ folder is already built and available, skip building it. 
    if [ -d "/ref/${params.genome}/STAR_index" ]; then
        ln -s ${params.inputdir}/../ref/${params.genome}/STAR_index STAR_index

    elif [ -d "${params.star_index}/STAR_index" ]; then
        ln -s ${params.star_index}/STAR_index STAR_index

    elif [ -w "/ref" ] && [ -w "/ref/${params.genome}" ]; then
        # If /ref and /ref/${params.genome} folders are writable, build STAR_index and copy there
        star_index="${params.inputdir}/../ref/${params.genome}/STAR_index"

    else 
        star_index="${params.star_index}/STAR_index"
    fi

    echo \${star_index}

 
    if [ ! -L STAR_index ]; then 
        source activate rna_env
        STAR --runMode genomeGenerate \
            --genomeDir \${star_index} \
            --runThreadN ${task.cpus} \
            --sjdbGTFfile ${gtf} \
            --genomeFastaFiles ${fasta} \
            --sjdbOverhang 100 
        conda deactivate
        ln -s \${star_index} STAR_index
    fi

    source activate samtools_env
    fasta="${fasta}"
    filename_no_ext=\${fasta%.fa}
    filename_no_ext=\${filename_no_ext%.fasta}
    if [ ! -f \${filename_no_ext}.chrom.sizes ]; then
        samtools faidx ${fasta}
        cut -f1,2 ${fasta}.fai >  \${filename_no_ext}.chrom.sizes
    fi
    ln -s \${filename_no_ext}.chrom.sizes genome.chrom.sizes
    conda deactivate
    """ 
}


