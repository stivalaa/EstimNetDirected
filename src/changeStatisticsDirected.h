#ifndef CHANGESTATISTICSDIRECTED_H
#define CHANGESTATISTICSDIRECTED_H
/*****************************************************************************
 * 
 * File:    changeStatisticsDirected.h
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: October 2017
 *
 * Functions to compute directed graph change statistics. Each
 * function takes a pointer to a digraph struct, and two node numbers
 * i and j and returns the value of the change statistic for adding
 * the arc i -> j.
 *
 * For change statistics dependent on a nodal attribute, there is
 * an additional parameter a which is the index of the attribute
 * to use.
 *
 * These functions are adapted from the original PNet code by Peng Wang:
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
 ****************************************************************************/

#include "utils.h"
#include "digraph.h"


/* typedef for change statistics function  */
typedef double (change_stats_func_t)(const digraph_t *g, uint_t i, uint_t j);

/* version for change statistics with nodal attribute */
typedef double (attr_change_stats_func_t)(const digraph_t *g, uint_t i, uint_t j, uint_t a);

/* version for change statistics with dyadic covariate */
/* for the moment just hte same as change_stats_func_t as treated specially,
   only used for GeoDistance for now */
typedef double (dyadic_change_stats_func_t)(const digraph_t *g, uint_t i, uint_t );

/* change statistics with pairs of nodal attributes (attribute interactions) */
typedef double (attr_interaction_change_stats_func_t)(const digraph_t *g, uint_t i, uint_t j, uint_t a, uint_t b);

/************************* Structural ****************************************/

double changeArc(const digraph_t *g, uint_t i, uint_t j);
double changeReciprocity(const digraph_t *g, uint_t i, uint_t j);
double changeSink(const digraph_t *g, uint_t i, uint_t j);
double changeSource(const digraph_t *g, uint_t i, uint_t j);
double changeInTwoStars(const digraph_t *g, uint_t i, uint_t j);
double changeOutTwoStars(const digraph_t *g, uint_t i, uint_t j);
double changeIsolates(const digraph_t *g, uint_t i, uint_t j);
double changeTwoPath(const digraph_t *g, uint_t i, uint_t j);
double changeTransitiveTriad(const digraph_t *g, uint_t i, uint_t j);
double changeCyclicTriad(const digraph_t *g, uint_t i, uint_t j);
double changeAltInStars(const digraph_t *g, uint_t i, uint_t j);
double changeAltOutStars(const digraph_t *g, uint_t i, uint_t j);
double changeAltKTrianglesT(const digraph_t *g, uint_t i, uint_t j);
double changeAltKTrianglesC(const digraph_t *g, uint_t i, uint_t j);
double changeAltKTrianglesD(const digraph_t *g, uint_t i, uint_t j);
double changeAltKTrianglesU(const digraph_t *g, uint_t i, uint_t j);
double changeAltTwoPathsT(const digraph_t *g, uint_t i, uint_t j);
double changeAltTwoPathsD(const digraph_t *g, uint_t i, uint_t j);
double changeAltTwoPathsU(const digraph_t *g, uint_t i, uint_t j);
double changeAltTwoPathsTD(const digraph_t *g, uint_t i, uint_t j);

/************************* Actor attribute (binary) **************************/

double changeSender(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeReceiver(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeInteraction(const digraph_t *g, uint_t i, uint_t j, uint_t a);

/********************* Actor attribute (categorical) *************************/

double changeMatching(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeMatchingReciprocity(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeMismatching(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeMismatchingReciprocity(const digraph_t *g, uint_t i, uint_t j, uint_t a);

/********************* Actor attribute (continuous) *************************/

double changeContinuousSender(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeContinuousReceiver(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeDiff(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeDiffReciprocity(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeDiffSign(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeDiffDirSR(const digraph_t *g, uint_t i, uint_t j, uint_t a);
double changeDiffDirRS(const digraph_t *g, uint_t i, uint_t j, uint_t a);


/********************* Actor attribute (set of categorical) *******************/

double changeJaccardSimilarity(const digraph_t *g, uint_t i, uint_t j, uint_t a);


/********************* Dyadic covariate (continuous) *************************/

double changeGeoDistance(const digraph_t *g, uint_t i, uint_t j);
double changeLogGeoDistance(const digraph_t *g, uint_t i, uint_t j);
double changeEuclideanDistance(const digraph_t *g, uint_t i, uint_t j);


/************ Actor attribute interaction (categorical) *********************/

double changeMatchingInteraction(const digraph_t *g, uint_t i, uint_t j,
                                 uint_t a, uint_t b);


/*************************** Other functions *********************************/

double calcChangeStats(const digraph_t *g, uint_t i, uint_t j,
                       uint_t n, uint_t n_attr, uint_t n_dyadic,
                       uint_t n_attr_interaction,
                       change_stats_func_t *change_stats_funcs[],
                       attr_change_stats_func_t *attr_change_stats_funcs[],
                       dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                       attr_interaction_change_stats_func_t 
                                        *attr_interaction_change_stats_funcs[],
                       uint_t attr_indices[],
                       uint_pair_t attr_interaction_pair_indices[],
                       const double theta[],
                       bool isDelete,
                       double changestats[]);


double jaccard_index(set_elem_e a[], set_elem_e b[], uint_t n);
double boundary_crossing_ratio(const digraph_t *g, uint_t i, uint_t a);

#endif /* CHANGESTATISTICSDIRECTED_H */

