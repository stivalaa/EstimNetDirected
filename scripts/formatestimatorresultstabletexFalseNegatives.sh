#!/bin/sh
#
# File:    formatestimatorresultstabletexFalseNegatives.sh
# Author:  Alex Stivala
# Created: November 2013
#
#
# sort rows of table output from makeestimatorresultstable.R and add
# LaTeX header/footer
# 
# Input is stdin
#
# Usage: formatestimatorresultstabletexFalseNegatives.sh estimatorresutlstablefile.txt
#
# E.g.:
#    formatestimatorresultstabletexFalseNegatives.sh  estimator_error_statistics.txt
#
# Output is to stdout
#
# Uses various GNU utils options on  echo &  etc.



# write_header() - write LaTeX table header for EstimNetDirected

write_header() {
  cat <<EOF
\begin{tabular}{rllrrrrrrrrr}
\hline
N &  Attributes &  Effect &  Bias &  RMSE &  \multicolumn{3}{c}{False negative rate (\%)}  & in C.I.    & Total     & Mean       & Total\\\\
  &             &         &       &       &   Estim. & \multicolumn{2}{c}{95\% C.I.}       &  (\%)      & networks  & runs       & runs per\\\\
  &             &         &       &       &   & lower & upper                              &            & converged & converged  & network\\\\
\hline
EOF
}

if [ $# -ne 1 ]; then
  echo "Usage: $0 error_statistics.txt" >&2
  exit 1
fi
infile=$1

write_header

grep EstimNetDirected $infile | sort -t\& -k4,4 -k3,3 -k1,1n -k2,2r  | awk -F\& -vOFS=\& '{printf("%s & %s & %s & %s & %5.4f & %5.4f & %2.0f & %2.0f & %2.0f & %2.0f & %d & %2.2f & %g\\\\\n", $1,$2,$3,$5,$6,$7,$8,$12,$13,$17,$15,$11,$16)}' | sed 's/ EstimNetDirected  &//g'

cat <<EOF 
\hline
\end{tabular}
EOF
