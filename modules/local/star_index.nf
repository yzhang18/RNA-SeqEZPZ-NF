/*
 * Define the "index" process
 * using tools STAR 
 */
process STAR_INDEX {
    label "star"
    publishDir params.star_index, mode: "copy", pattern: "STAR_index"

    input:
    val (fasta) 
    val (gtf) 

    output:
    path ("*.chrom.sizes"), emit: chr
    path ("STAR_index"), emit: index

    script:
    """ 
    source activate rna_env
    STAR --runMode genomeGenerate \
        --genomeDir STAR_index \
        --runThreadN ${task.cpus} \
        --sjdbGTFfile ${gtf} \
        --genomeFastaFiles ${fasta} \
        --sjdbOverhang 100 
    conda deactivate

    source activate samtools_env
    fasta="${fasta}"
    filename_no_ext=\${fasta%.fa}
    filename_no_ext=\${filename_no_ext%.fasta}
    if [ ! -f \${filename_no_ext}.chrom.sizes ]; then
        samtools faidx ${fasta}
    fi
    cut -f1,2 ${fasta}.fai >  \${filename_no_ext}.chrom.sizes
    cut -f1,2 ${fasta}.fai > genome.chrom.sizes
    conda deactivate
    """ 
}


