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


/************************* Structural ****************************************/

double changeArc(const digraph_t *g, uint_t i, uint_t j);
double changeReciprocity(const digraph_t *g, uint_t i, uint_t j);
double changeAltInStars(const digraph_t *g, uint_t i, uint_t j);
double changeAltOutStars(const digraph_t *g, uint_t i, uint_t j);
double changeAltKTrianglesT(const digraph_t *g, uint_t i, uint_t j);
double changeAltKTrianglesC(const digraph_t *g, uint_t i, uint_t j);
double changeAltTwoPathsT(const digraph_t *g, uint_t i, uint_t j);
double changeAltTwoPathsD(const digraph_t *g, uint_t i, uint_t j);
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

/********************* Dyadic covariate (continuous) *************************/

double changeGeoDistance(const digraph_t *g, uint_t i, uint_t j);
double changeLogGeoDistance(const digraph_t *g, uint_t i, uint_t j);

#endif /* CHANGESTATISTICSDIRECTED_H */

