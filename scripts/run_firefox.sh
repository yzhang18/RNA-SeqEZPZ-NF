#!/bin/bash
set -x
singularity exec $img_dir/$img_name bash -c "source activate firefox_env; \
	firefox --no-remote --new-window -P \"default\" http://127.0.0.1:$port_num"
