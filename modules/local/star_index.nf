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
        --sjdbGTFfile /ref/${params.genome}/${gtf} \
        --genomeFastaFiles /ref/${params.genome}/${fasta} \
        --sjdbOverhang 100 
    conda deactivate

    source activate samtools_env
    if [ ! -f /ref/${params.genome}/${fasta} ]; then
        samtools faidx /ref/${params.genome}/${fasta}
    fi
    cut -f1,2 /ref/${params.genome}/${fasta}.fai >  ${fasta}.chrom.sizes
    conda deactivate
    """ 
}


