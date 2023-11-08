run on change_input. run with hg38 that already has star_index in image directory done.

run on change_input. run with ref_fa and ref_gtf path to danRer11 error in run_sartools I think due to human samples

run on change_input. run with ref_fa and ref_gtf path to a673_hap1_0

rerun above after adding khmer to calculate genome

rerun above after fixing fail in calculating readlength and not passing readlength to star_pass2, also fixed failing because there is no chrom sizes but star index exists. update actually it is not failing because there is no chrom size but now I'm checking for chrom sizes and fasta index before skipping running star_index_simg

run on change_input. with genome=danRer11 and ref_fa and ref_gtf
