/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NCH RNASEQ NEXTFLOW CONFIG FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT INFORMATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.inputdir = "project_ex"
params.input = "$params.inputdir/samples.txt"
//params.reads = "$params.inputdir/fastq/*_R{1,2}_*.fastq.gz"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUT DIRECTORIES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.outdir = "$params.inputdir"
params.star_index = "$params.outdir/outputs"
params.merged_fastq = "$params.outdir/outputs/merged_fastq"
params.trim = "$params.outdir/outputs/trim"
params.fastqc = "$params.outdir/outputs/fastqc_rslt"
params.star_pass1 = "$params.outdir/outputs/STAR_2pass/Pass1"
params.star_pass2 = "$params.outdir/outputs/STAR_2pass/Pass2"
params.star_genome_for_pass2 = "$params.outdir/outputs/STAR_2pass/GenomeForPass2"
params.bw_files = "$params.outdir/outputs/bw_files"
params.sampledir = "$params.outdir/outputs/samples"
params.multiqc = "$params.outdir/outputs"
params.shiny_app = "$params.outdir/outputs/shiny_app"
params.logdir = "$params.outdir/outputs/logs"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.ref = "$params.inputdir/ref"
params.genome = "hg19"

if ( "$params.genome" == "hg19" )
{
    params.fasta = "human_g1k_v37_decoy.fasta"
}
else if ( "$params.genome" == "hg38" )
{
    params.fasta = "hg38.analysisSet.fa" 
}
else
{
    params.fasta = "*.fasta"
}
params.gtf = "*.gtf"



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OTHER PARAMETERs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.padj = "0.05"
params.prjdir = ""
params.email = ""
params.batch_adjust = "yes"
params.ncpus_trim = 4
params.ncpus_star = 20 

params.debug = 'no'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
WorkflowMain.initialise(workflow, params, log)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN RNASEQ PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RNASEQ } from './workflows/run_rnaseq_full'

workflow {
    RNASEQ ()
}

