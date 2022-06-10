/*****************************************************************************
 * 
 * File:    testRandomIntRange.c
 * Author:  Alex Stivala
 * Created: June 2022
 *
 * Generate random bounded integers to test for uniformity etc.
 *
 *
 * Usage:  testRandomIntRange <n> <num>
 *
 * Generate <num> random integers in range 0...<n>-1
 *
 * Output is to stdout, one integer per line
 * 
 *
 ****************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "utils.h"


int main(int argc, char *argv[]) 
{
  uint_t n, num, i, r;

  init_prng(0); /* initialize pseudorandom number generator */

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <n> <num>\n", argv[0]);
    exit(1);
  }
  n = atol(argv[1]);
  num = atol(argv[2]);

  for (i = 0; i < num; i++) {
    r = int_urand(n);
    assert(r >= 0 && r < n);
    printf("%u\n", r);
  }
  exit(0);
}

