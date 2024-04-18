#!/bin/bash
set -x
singularity exec \
	--bind $proj_dir:/mnt \
	--bind $img_dir/scripts:/scripts -B ~/.Xauthority \
	$bind_filepath \
	--bind $proj_dir/mypipe:/hostpipe \
        $img_dir/$img_name /bin/sh /scripts/app_simg.sbatch &> run_shiny_analysis.out &

# Check to make sure shiny app is already at listening point
# before continuing
# Put this in because sometimes it takes time to load libraries
while true; do 
        # get the last line of run_shiny_analysis.out
        last_line=$(tail -n 1 $proj_dir/run_shiny_analysis.out)
        # exit if shiny fail to load
        if grep -q "Execution halted" "run_shiny_analysis.out"; then exit 1; fi
        # check if it contains listening
        if [[ $last_line == *"Listening"* ]]; then
                break
        fi
	sleep 5
done

# listen to shiny app
cd $proj_dir
eval "$(cat mypipe)"
