#!/bin/sh
##
## Run the R scripts to plot fourcycle dependence figures.
##
## The CYPATH directory containing cypath executable and Perl scripts
## (transgrh.pl etc.) must be in the PATH as this R script uses
## system2() to call them in order to find (not just count)
## four-cycles. Have to edit the transgrh.pl script to add a
## semicolon on the end of line 38 as it does not work otherwise.
##
## For CYPATH see:
##
##    http://research.nii.ac.jp/~uno/code/cypath.html
##    http://research.nii.ac.jp/~uno/code/cypath11.zip
## 
##    Uno, T., & Satoh, H. (2014, October). An efficient algorithm for
##    enumerating chordless cycles and chordless paths. In International
##    Conference on Discovery Science (pp. 313-324). Springer, Cham.
##
##

CYPATH_DIR=${DOCUMENTS}/USI/CYPATH
export PATH=${PATH}:${CYPATH_DIR}
Rscript plot_lattice_bipartite_fourcycles_dependence_figures.R
Rscript plot_example_fourcycles_dependence_figures.R
