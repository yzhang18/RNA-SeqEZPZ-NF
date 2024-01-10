#!/bin/bash

source activate sartools_env;
set -x
cd /scripts
#hostfilepath="\"$hostfilepath\""
hostfilepath="$hostfilepath"
R -e "hostfilepath='$hostfilepath'; max.nsamples=$max_nsamples; shiny::runApp(port=$port_num)"
conda deactivate
