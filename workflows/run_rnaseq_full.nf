/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters
 if (params.inputdir) { params.input = "$params.inputdir/samples.txt" } else { exit 1, 'Input directory not specified!' }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULE: Loaded from modules/local/
include { MERGE_FASTQ   } from '../modules/local/merge_fastq'
include { TRIM_FASTQC   } from '../modules/local/trim_fastqc'
include { MULTIQC1      } from '../modules/local/multiqc1'
include { MULTIQC2      } from '../modules/local/multiqc2'
include { STAR_INDEX    } from '../modules/local/star_index'
include { STAR_PASS1    } from '../modules/local/star_pass1'
include { GENOME_SIZE   } from '../modules/local/genome_size'
include { STAR_PASS2    } from '../modules/local/star_pass2'
include { COMBINEBW     } from '../modules/local/combinebw'
include { FEATURE_COUNTS} from '../modules/local/feature_counts'
include { SARTOOLS      } from '../modules/local/sartools'
include { GET_SAMPLE    } from '../modules/local/get_sample'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// SUBWORKFLOW:
include { INPUT_CHECK    } from '../subworkflows/nf-core/input_check'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RNASEQ {
    ch_samples = file(params.input) 
    ch_fasta = Channel.value (params.fasta)
    ch_gtf = Channel.value(params.gtf)
    // ch_fasta.view()

    /*
     * MODULE: STAR_INDEX
     */
    STAR_INDEX (
        ch_fasta,
        ch_gtf
    )
    ch_chr = STAR_INDEX.out.chr
    ch_index = STAR_INDEX.out.index


    /*
     * MODULE: MERGE_FASTQ
     */
    MERGE_FASTQ (
        ch_samples
    )

    ch_merged_samples = MERGE_FASTQ.out.ch_merged_samples

    /*
     * SUBWORKFLOW: INPUT_CHECK
     * 
     * Read in samples.csv, validate and stage input files
     */
    INPUT_CHECK (
        ch_merged_samples
    )

    ch_fastq = INPUT_CHECK.out.reads


    /*
     * MODULE: trim the reads and run fastqc
     */
    TRIM_FASTQC (
        ch_fastq,
        ch_chr
    )
    
    reads = TRIM_FASTQC.out.reads
    report = TRIM_FASTQC.out.report
    html = TRIM_FASTQC.out.html
    zip = TRIM_FASTQC.out.zip
    ch_log = TRIM_FASTQC.out.log

    /*
     * MODULE: run multiqc for the first time 
     */
    MULTIQC1(
        ch_log.collect()
    )    
    readlength = MULTIQC1.out.readlength
 
    /*
     * MODULE: STAR_PASS1 
     */
    STAR_PASS1 ( reads, ch_index )

    sjFiles = STAR_PASS1.out.sjdb.collect()


    /*
     * MODULE: calculate genome size
     */
    GENOME_SIZE (ch_fasta, readlength)
    size = GENOME_SIZE.out.size

    /*
     * MODULE: STAR_PASS2
     */
    STAR_PASS2 ( reads, sjFiles, size, ch_index)

    ch_bam_sorted = STAR_PASS2.out.bam_sorted
    ch_bw = STAR_PASS2.out.bw
    .map {
        meta, fastq ->
            [ meta.groupid, fastq ] }
    .groupTuple(by: [0])


    /*
     * MODULE: COMBINEBW
     */
    COMBINEBW ( ch_chr, ch_bw)

    /*
     * MODULE: FEATURE_COUNT
     */
    FEATURE_COUNTS (ch_bam_sorted, ch_gtf)

    ch_counts = FEATURE_COUNTS.out.counts
    .map {
        meta, count_file ->
            count_file  }
    .collect()

    /*
     * MODULE: SARTOOLS
     */
    SARTOOLS (ch_counts)
    ch_log = SARTOOLS.out.log

    /*
     * MODULE: MULTIQC2
     * rerun MULTIQC at the end of the pipeline
     */
    MULTIQC2 (ch_log)

}


