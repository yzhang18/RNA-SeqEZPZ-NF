#!/bin/bash

source activate sartools_env;
set -x
cd /scripts
R -e "shiny::runApp(port=$port_num)"
conda deactivate