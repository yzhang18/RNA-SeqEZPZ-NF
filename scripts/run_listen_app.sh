#!/bin/bash
set -x
set -v

# listen to shiny app
cd $proj_dir
eval "$(cat mypipe)" &

SINGULARITYENV_port_num=$port_num \
        SINGULARITYENV_hostfilepath=$filepath \
        SINGULARITYENV_max_nsamples=$max_nsamples \
        SINGULARITYENV_img_dir=$img_dir \
	singularity exec \
		--bind $proj_dir:/mnt \
		--bind $img_dir/scripts:/scripts -B ~/.Xauthority \
		$bind_filepath \
		--bind $proj_dir/mypipe:/hostpipe \
		--bind $img_dir/ref:/ref \
		--bind $img_dir:/img_dir \
        	$img_dir/$img_name /bin/sh /scripts/app_simg.sbatch

