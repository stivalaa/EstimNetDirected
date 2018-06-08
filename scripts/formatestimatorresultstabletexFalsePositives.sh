#!/bin/sh
#
# File:    formatestimatorresultstabletexFalsePositives.sh
# Author:  Alex Stivala
# Created: December 2013
#
#
# sort rows of table output from makeMLEresultstable.R and add
# LaTeX header/footer
# 
# Input is stdin
#
# Usage: formatestimatorresultstabletexFalsePositvies.sh estimatorresutlstablefile.txt
#
# E.g.:
#    formatestimatorresultstabletex.sh  estimator_error_false_positives_statistics.txt
#
# Output is to stdout
#
# Uses various GNU utils options on  echo &  etc.


# write_header() - write LaTeX table header for EstimNetDirected
write_header() {
  cat <<EOF 
\begin{tabular}{rllrrrrrrrrr}
\hline
N &  Attributes &  Effect &  Bias &  RMSE &  \multicolumn{3}{c}{False positive rate (\%)}  & in C.I.    & Total     & Mean       & Total\\\\
  &             &         &       &       &   Estim. & \multicolumn{2}{c}{95\% C.I.}       &  (\%)      & networks  & runs       & runs per\\\\
  &             &         &       &       &   & lower & upper                              &            &           & converged  & network\\\\
\hline
EOF
}


if [ $# -ne 1 ]; then
  echo "Usage: $0 error_statistics.txt" >&2
  exit 1
fi
infile=$1



write_header

# the --posix flag on awk (gawk) stop it converting NaN and Inf to 0
# then sed 's/nan/--/g' actuall converts to --- as it has written nan as -nan
grep EstimNetDirected $infile | sort -t\& -k4,4 -k3,3 -k1,1n -k2,2r  | awk --posix -F\& -vOFS=\& '{ if ($9 != " NA ") printf("%s & %s & %s & %s & %5.4f & %5.4f & %2.0f & %2.0f & %2.0f & %2.0f & %d & %2.2f & %g\\\\\n",$1,$2,$3,$5,$7,$8,$9,$13,$14,$18,$16,$12,$17)}' | sed 's/nan/--/g'  | sed 's/ EstimNetDirected  &//g'  

cat <<EOF
\hline
\end{tabular}
EOF
