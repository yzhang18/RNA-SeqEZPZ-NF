#!/bin/bash
echo $run
if [[ $run == "debug" ]];then
	set -x
fi
singularity exec $img_dir/$img_name bash -c "source activate firefox_env; \
	firefox --no-remote --new-window -CreateProfile shiny_${port_num}; \
	firefox --no-remote --new-window -P shiny_${port_num} http://127.0.0.1:$port_num"

# remove the shiny profile in firefox
cd ~/.mozilla/firefox
ls | grep shiny_${port_num} | xargs rm -r
