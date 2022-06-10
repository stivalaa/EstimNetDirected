#!/usr/bin/env Rscript
##############################################################################
#
# File:    checkUniform.R
# Author:  Alex Stivala
# Created: June 2022
# checkUniform.py - use chi-squared test to check for uniform distribution
#
# Reads output of testRandom on stdin and does chi-squared test for
# discrete uniform distribution.
#
# (Note that even using the basic rand() % n which is biased, this test
# will likely pass, the bias is just not large enough to detect this way,
# but at least this way we are doing some testing, including in the C code
# to make sure the range is not actually wrong (assert) ec.)
#
#
# Usage: Rscript checkUniform.R n
#  where n is the number of possible integers i.e. range is 0 .. n-1
#
##############################################################################

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Usage: Rscript checkUniform.R n\n")
  quit(save="no", status=2)
}
n <- as.integer(args[1])

con <- file("stdin")
x <- scan(con)
# https://stackoverflow.com/questions/1617061/include-levels-of-zero-count-in-result-of-table
frequencies <- table(factor(x, levels = 0 : (n-1)))
##print(frequencies)#XXX
testresult <- chisq.test(frequencies)
print(testresult)
if (testresult$p.value < 0.05) {
  print("Rejected uniform distribution")
  quit(status = 1)
} else {
  print("Consistent with uniform distribution")
}
close(con) # stop warning 1: In q : closing unused connection 3 (stdin)
quit(status = 0)

