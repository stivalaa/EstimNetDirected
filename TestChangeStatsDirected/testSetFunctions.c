/*****************************************************************************
 * 
 * File:    testSetFunctions.c
 * Author:  Alex Stivala
 * Created: October 2017
 *
 * Test set (of categorical) parsing and Jaccard similarity etc.
 *
 *
 * Usage:  testSetFunctions  <infilename>
 *
 * Reads single set attribute values
 * from <infilename> in EstimNetDirected format.
 *
 * Output is to stdout, each line is:
 *   i j value
 * where a and j are indices (rows in input file, with 0 as first row after
 * header line) and value is Jaccard index for the two sets at rows i and j
 * in the input.
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
#include "digraph.h"
#include "changeStatisticsDirected.h"

#define MAX_VALS  2000
#define NUM_TESTS 1000

int main(int argc, char *argv[]) 
{
  char buf[1024];
  uint_t i,j,k;
  char *set_filename = NULL;
  FILE *file           = NULL;
  int pass;
  uint_t size = 0;
  set_elem_e *set[MAX_VALS];
  uint_t numvals = 0;
  double sim;

  srand(time(NULL));

  if (argc != 2) {
    fprintf(stderr, "Usage: %s <infilename>\n", argv[0]);
    exit(1);
  }
  set_filename = argv[1];

  for (pass = 0; pass < 2; pass++) {
    numvals = 0;
    if (!(file = fopen(set_filename, "r"))) {
      fprintf(stderr, "error opening file %s (%s)\n", 
              set_filename, strerror(errno));
      return -1;
    }
    if (!fgets(buf, sizeof(buf)-1, file)) {
      fprintf(stderr, "ERROR: could not read header line in set attributes file %s (%s)\n",
              set_filename, strerror(errno));
      return -1;
    }
    if (!fgets(buf, sizeof(buf)-1, file)) {
      if (!feof(file)) {
        fprintf(stderr, "ERROR: attempting to read first line of set attributes in file %s (%s)\n",
                set_filename, strerror(errno));
        return -1;
      }
    }
    while (!feof(file)) {
      if (numvals >= MAX_VALS) {
        fprintf(stderr, "input is too large, limited to %d\n", MAX_VALS);
        exit(1);
      }
      rstrip(buf);
      if (parse_category_set(buf, pass == 0, &size, set[numvals]) < 0) {
        fprintf(stderr, "ERROR parsing set\n");
        exit(1);
      }
      numvals++;
      
      if (!fgets(buf, sizeof(buf)-1, file)) {
        if (!feof(file)) {
          fprintf(stderr, "ERROR: attempting to read set attributes in file %s (%s)\n",
                  set_filename, strerror(errno));
          return -1;
        }
      }
    }
    fclose(file);
    if (pass == 0) {
      for (k = 0; k < numvals; k++) {
        set[k] = (set_elem_e *)safe_malloc(size * sizeof(set_elem_e));
      }
    }
  }

  for (k = 0; k < NUM_TESTS; k++) {
    i = rand() % numvals;
    j = rand() % numvals;
    sim =  jaccard_index(set[i], set[j], size);
    assert(sim >= 0 && sim <= 1);
    printf("%u %u %f\n", i , j, sim);
  }
  exit(0);
}
  
