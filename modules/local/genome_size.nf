/*
 * define genome_size process
 *
 */
process GENOME_SIZE {
    label "hi_mem_cpus"
    input:
    val (fasta)
    path (readlength)

    output:
    path ("gsize.txt"), emit: size

    script:
    """
    # calculating effective genome size
    source activate khmer_env
    set -x

    read readlength < ${readlength}
    unique-kmers.py -k \$readlength /ref/${params.genome}/${fasta} &> ./kmers.txt
    genome_size=\$(tail -n1 ./kmers.txt | cut -f2 -d:)
    echo \$genome_size > gsize.txt
    """
}


