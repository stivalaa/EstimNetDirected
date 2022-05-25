#ifndef CHANGESTATISTICSUNDIRECTED_H
#define CHANGESTATISTICSUNDIRECTED_H
/*****************************************************************************
 * 
 * File:    changeStatisticsUndirected.h
 * Author:  Alex Stivala
 * Created: January 2022
 *
 * Functions to compute graph change statistics for undirected graphs. Each
 * function takes a pointer to a graph struct, and two node numbers
 * i and j and returns the value of the change statistic for adding
 * the edge i -- j.
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

#include "utils.h"
#include "graph.h"
#include "changeStatisticsTypes.h"

/************************* Utility functions *********************************/
ulong_t num_s_stars(const graph_t *g, uint_t v, ulong_t s);
ulong_t change_s_stars(const graph_t *g, uint_t v, ulong_t s);

  
/************************* Structural ****************************************/

double changeEdge(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltStars(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltTwoPaths(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltKTriangles(graph_t *g, uint_t i, uint_t j, double lambda);
double changeFourCycles(graph_t *g, uint_t i, uint_t j, double lambda);
double changeThreePaths(graph_t *g, uint_t i, uint_t j, double lambda);

/************************* Actor attribute (binary) **************************/

double changeActivity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);


#endif /* CHANGESTATISTICSUNDIRECTED_H */

