/*****************************************************************************
 * 
 * File:    changeStatisticsGeneral.c
 * Author:  Alex Stivala
 * Created: January 2022
 *
 * Functions to compute graph change statistics, that apply to both
 * directed and undirected graphs. Each
 * function takes a pointer to a digraph struct, and two node numbers
 * i and j and returns the value of the change statistic for adding
 * the edge i -- j or arc i -> j or j <- i; the functions in this 
 * module are for those where these cases are all the same.
 *
 * Also takes lambda (decay) parameter which is only used for
 * some statistics ("alternating" statistics).
 *
 * For change statistics dependent on a nodal attribute, there is
 * an additional parameter a which is the index of the attribute
 * to use.
 *
 * On some functions there is additionally a parameter indicating when
 * the change statistic is being computed as part of a delete (rather
 * than add) move, which can be used for some implementations that can
 * be more easily implemented with this information. However in
 * general it is more elegant and simpler to compute the statistic for
 * adding the arc (for delete moves the value returned is just
 * negated, and the change statistic function does not depend on or
 * need to use this information at all).
 *
 * Some of these functions are adapted from the original PNet code by Peng Wang:
 *
 *   Wang P, Robins G, Pattison P. PNet: A program for the simulation and
 *   estimation of exponential random graph models. University of
 *   Melbourne. 2006.
 *
 * And for the definitions of the change statistics:
 * 
 *   Robins, G., Pattison, P., & Wang, P. (2009). Closure, connectivity and
 *   degree distributions: Exponential random graph (p*) models for
 *   directed social networks. Social Networks, 31(2), 105-117.
 * 
 *   Snijders, T. A., Pattison, P. E., Robins, G. L., & Handcock,
 *   M. S. (2006). New specifications for exponential random graph
 *   models. Sociological methodology, 36(1), 99-153.
 * 
 * And also generally:
 * 
 *   Lusher, D., Koskinen, J., & Robins, G. (Eds.). (2013). Exponential
 *   random graph models for social networks: Theory, methods, and
 *   applications. New York, NY: Cambridge University Press.
 * 
 * especially Ch. 6:
 *
 *   Koskinen, J., & Daraganova, G. (2013). Exponential random graph model
 *   fundamentals. In "Exponential random graph models for social networks:
 *   Theory, methods, and applications." (pp. 49-76). New York, NY:
 *   Cambridge University Press.
 *
 * As well as the statnet ergm terms, and references for specific
 * change statistics where indicated.
 *
 *
 * Do NOT compile with -ffast-math on gcc as we depend on IEEE handling of NaN
 *
 ****************************************************************************/

#include <math.h>
#include <assert.h>
#include "changeStatisticsGeneral.h"


/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/

/*
 * signum function, returns -1 for negative x, +1 for positive x, else 0
 */
static double signum(double x)
{
  /* https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c */
  return (0 < x) - (x < 0);
}


/*
 * Size of the intersection of two sets.  Each set is represented by
 * array of set_elem_e indicating NA, ABSENT or PRESENT for each
 * element, size of intersection is number of array indices where both
 * are PRESENT. Note we ignore NA here.
 *
 * Parameters:
 *      a - set
 *      b - set
 *      n - size of arrays a and b
 *
 * Return value:
 *      number of elements that are in both sets a and b.
 */
static uint_t set_intersection_size(set_elem_e a[], set_elem_e b[], uint_t n)
{
  uint_t i;
  uint_t count = 0;

  for (i = 0; i < n; i++) {
    if (a[i] == SET_ELEM_PRESENT && b[i] == SET_ELEM_PRESENT) {
      count++;
    }
  }
  return count;
}

/*
 * Size of the union of two sets.  Each set is represented by
 * array of set_elem_e indicating NA, ABSENT or PRESENT for each
 * element, size of union is number of array indices where either (or both)
 * are PRESENT. Note we ignore NA here.
 *
 * Parameters:
 *      a - set
 *      b - set
 *      n - size of arrays a and b
 *
 * Return value:
 *      number of elements that are in either (or both) of the sets a and b.
 */
static uint_t set_union_size(set_elem_e a[], set_elem_e b[], uint_t n)
{
  uint_t i;
  uint_t count = 0;

  for (i = 0; i < n; i++) {
    if (a[i] == SET_ELEM_PRESENT || b[i] == SET_ELEM_PRESENT) {
      count++;
    }
  }
  return count;
}

/*
 * Jaccard index (similarity) for two sets. The Jaccard index is the
 * size of the intersection over the size of the union. (If a and b 
 * are both empty it is defined as 1).
 *
 * Parameters:
 *      a - set
 *      b - set
 *      n - size of arrays a and b
 *
 * Return value:
 *      Jaccard coefficient (similarity) of the two sets a and b
 *                                              
 * Has external linkage for use in unit tests.
 */
double jaccard_index(set_elem_e a[], set_elem_e b[], uint_t n)
{
  uint_t intersection_size, union_size;

  intersection_size = set_intersection_size(a, b, n);
  union_size = set_union_size(a, b, n);
  return (union_size == 0) ? 1 : (double)intersection_size / (double)union_size;
}


/*****************************************************************************
 *
 * change statistics functions
 *
 ****************************************************************************/


/********************* Actor attribute (continuous) *************************/

/*
 * Change statistic for continuous diff sign (sign of difference of attribute)
 * for attr_i - attr_j (so +1 when sending node has higher attribute value and -1
 * when receiving node has higher attribute value).
 */
double changeDiffSign(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /* unused parameter */
  if (isnan(g->contattr[a][i]) || isnan(g->contattr[a][j]))
    return 0;
  else
    return signum(g->contattr[a][i] - g->contattr[a][j]);
}


/***************** Actor attribute (set of categorical) ********************/


/*
 * Change statistic for set Jaccard similarity
 */
double changeJaccardSimilarity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete)
{
  (void)isDelete; /* unused parameter */    
  /* For NA values all elements of set are set to NA so just check first */
  if (g->setattr[a][i][0] == SET_ELEM_NA || g->setattr[a][j][0] == SET_ELEM_NA)
    return 0;
  else
    return jaccard_index(g->setattr[a][i], g->setattr[a][j],
                         g->setattr_lengths[a]);
}

