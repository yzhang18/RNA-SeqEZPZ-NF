#!/bin/bash
set -x
# deleting existing profile to prevent it from being locked.
rm -rf ~/.config/google-chrome/Singleton*
singularity exec $img_dir/$img_name bash -c "source activate chrome_env; \
	google-chrome http://127.0.0.1:$port_num"
