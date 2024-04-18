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
        	$img_dir/$img_name /bin/sh /scripts/app_simg.sbatch
echo "Right after shiny app" > test.txt
# Check to make sure shiny app is already at listening point
# before continuing
# Put this in because sometimes it takes time to load libraries
#while true; do 
#        # get the last line of run_shiny_analysis.out
#        last_line=$(tail -n 1 $proj_dir/run_shiny_analysis.out)
#        # exit if shiny fail to load
#        if grep -q "Execution halted" "run_shiny_analysis.out"; then exit 1; fi
#        # check if it contains listening
#        if [[ $last_line == *"Listening"* ]]; then
#                break
#        fi
#	sleep 5
#done

