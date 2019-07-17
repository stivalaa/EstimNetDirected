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
 ****************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include "digraph.h"
#include "changeStatisticsDirected.h"


int main(int argc, char *argv[]) 
{
  char buf[1024];
  uint_t i,j;
  char *set_filename = NULL;
  FILE *file           = NULL;
  int pass;
  uint_t size = 0;
  set_elem_e *set = NULL;
 
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <infilename>\n", argv[0]);
    exit(1);
  }
  set_filename = argv[1];

  for (pass = 0; pass < 2; pass++) {
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
      rstrip(buf);
      if (parse_category_set(buf, pass == 0, &size, set) < 0) {
        fprintf(stderr, "ERROR parsing set\n");
        exit(1);
      }
      
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
      set = (set_elem_e *)safe_malloc(size * sizeof(set_elem_e));
    }
  }
  exit(0);
}
  
