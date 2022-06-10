##############################################################################
#
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
# File:    checkSetFunctions.py
# Author:  Alex Stivala
# Created: June 2022
#
##############################################################################

con <- file("stdin")
x <- scan(con)
testresult <- chisq.test(table(x))
print(testresult)
if (testresult$p.value < 0.05) {
  print("Rejected uniform distribution")
  quit(status = 1)
} else {
  print("Consistent with uniform distribution")
}
close(con) # stop warning 1: In q : closing unused connection 3 (stdin)
quit(status = 0)

