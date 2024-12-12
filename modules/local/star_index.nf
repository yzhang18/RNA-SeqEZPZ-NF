/*
 * Define the "index" process
 * using tools STAR 
 */
process STAR_INDEX {
    label "star"
    publishDir params.logdir, mode: "copy", pattern: "star_index.out"
    //publishDir params.star_index, mode: "copy", pattern: "STAR_index"

    input:
    val (fasta) 
    val (gtf) 

    output:
    path ("genome.chrom.sizes"), emit: chr
    path ("STAR_index"), emit: index
    path ("star_index.out"), emit: log

    script:
    """
    set -x 

    # Set the bash variable fasta.
    fasta=${fasta}

    # Get the folder where the fasta file is.
    fasta_dir=\$(dirname \${fasta})

    # The first location where STAR_index could be.
    star_index_1=\${fasta_dir}/${params.genome}/STAR_index

    # The second location where STAR_index could be
    star_index_2=${params.inputdir}/${params.genome}/STAR_index

    # The variable to hold the location of STAR_index if it's not built.
    star_index=

    # If STAR_index/ folder is already built and available, just copy it here.
    if [ -d \${star_index_1} ] && [ -z \$(ls -A \${star_index_1} >/dev/null 2>&1) ]; then
        ln -s \${star_index_1} STAR_index
    elif [ -d \${star_index_2} ]; then
        ln -s \${star_index_2} STAR_index
    fi

    # If STAR_index isn't built, build it in its location and then copy it here. 
    if [ ! -L STAR_index ]; then 
        if [ -w \${fasta_dir} ]; then
            star_index=\${star_index_1}
        else 
            star_index=\${star_index_2}
        fi  

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
    filename_no_ext=\${fasta%.fa}
    filename_no_ext=\${filename_no_ext%.fasta}
    if [ ! -f \${filename_no_ext}.chrom.sizes ]; then
        samtools faidx ${fasta}
        cut -f1,2 ${fasta}.fai >  \${filename_no_ext}.chrom.sizes
    fi
    ln -s \${filename_no_ext}.chrom.sizes genome.chrom.sizes
    conda deactivate
    
    cat .command.log > star_index.out
    """ 
}


