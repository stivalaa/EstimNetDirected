#!/bin/sh
#
# File:    formatestimatorresultstabletex.sh
# Author:  Alex Stivala
# Created: November 2013
#
#
# sort rows of table output from makeMLEresultstable.R and add
# LaTeX header/footer
# 
# Input is stdin
#
# Usage: formatestimatorresultstabletex.sh estimatorresutlstablefile.txt
#
# E.g.:
#    formatestimatorresultstabletex.sh  estimator_error_statistics.txt
#
# Output is to stdout.
#
# Uses various GNU utils options on  echo &  etc.


# write_header() - write LaTeX table header for EstimNetDirected table 
# to stdout
write_header() {
  cat <<EOF
\begin{tabular}{rllrrrrrrr}
\hline
N &  Attributes  &  Effect &  Bias &  RMSE &  Mean           & Std. dev.   & Total    &   Mean     &  Total \\\\
  &              &         &       &       &  standard       &   estimate  & networks & runs       &  runs per \\\\
  &              &         &       &       &  error          &             &          & converged  &  network  \\\\
\hline
EOF
}

if [ $# -ne 1 ]; then
  echo "Usage: $0 error_statistics.txt" >&2
  exit 1
fi
infile=$1


write_header

grep EstimNetDirected $infile | sort -t\& -k4,4 -k3,3 -k1,1n -k2,2r  | awk -F\& -vOFS=\& '{printf("%s & %s & %s & %s & %5.4f & %5.4f & %5.4f & %5.4f & %d & %0.2f & %g \\\\\n", $1,$2,$3,$5,$6,$7,$9,$10,$15,$11,$16)}' | sed 's/  EstimNetDirected  &//g'  

cat <<EOF
\hline
\end{tabular}
EOF
