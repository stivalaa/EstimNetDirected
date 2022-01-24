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
 ****************************************************************************/

#include "utils.h"
#include "graph.h"
#include "changeStatisticsTypes.h"


/************************* Structural ****************************************/

double changeArc(graph_t *g, uint_t i, uint_t j, double lambda);
double changeReciprocity(graph_t *g, uint_t i, uint_t j, double lambda);
double changeSink(graph_t *g, uint_t i, uint_t j, double lambda);
double changeSource(graph_t *g, uint_t i, uint_t j, double lambda);
double changeInTwoStars(graph_t *g, uint_t i, uint_t j, double lambda);
double changeOutTwoStars(graph_t *g, uint_t i, uint_t j, double lambda);
double changeTransitiveTriad(graph_t *g, uint_t i, uint_t j, double lambda);
double changeCyclicTriad(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltInStars(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltOutStars(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltKTrianglesT(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltKTrianglesC(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltKTrianglesD(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltKTrianglesU(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltTwoPathsT(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltTwoPathsD(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltTwoPathsU(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltTwoPathsTD(graph_t *g, uint_t i, uint_t j, double lambda);
double changeLoopInteraction(graph_t *g, uint_t i, uint_t j, double lambda);

/************************* Actor attribute (binary) **************************/

double changeSender(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeReceiver(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);


/********************* Actor attribute (categorical) *************************/

double changeMatchingReciprocity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeMismatchingReciprocity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeMismatchingTransitiveTriad(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeMismatchingTransitiveTies(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);

/********************* Actor attribute (continuous) *************************/

double changeContinuousSender(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeContinuousReceiver(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeDiffReciprocity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeDiffSign(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeDiffDirSR(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);
double changeDiffDirRS(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete);



#endif /* CHANGESTATISTICSDIRECTED_H */

