#!/bin/bash

source activate sartools_env;
cd /scripts
R -e "if(length(.libPaths())>1) .libPaths(.libPaths()[-1]); shiny::runApp(port=$port_num)"
conda deactivate
