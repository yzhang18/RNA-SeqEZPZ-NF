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
    chr_info=${fasta%.*}
    if [ ! -f ${chr_info}.chrom.sizes ]; then
        samtools faidx ${fasta}
    fi
    cut -f1,2 ${chr_info}.fai >  ${chr_info}.chrom.sizes
    conda deactivate
    """ 
}


