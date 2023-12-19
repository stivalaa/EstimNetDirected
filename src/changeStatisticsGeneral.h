#ifndef CHANGESTATISTICSGENERAL_H
#define CHANGESTATISTICSGENERAL_H
/*****************************************************************************
 * 
 * File:    changeStatisticsGeneral.h
 * Author:  Alex Stivala
 * Created: January 2022
 *
 * Functions to compute graph change statistics, that apply to both
 * directed and undirected graphs. Each
 * function takes a pointer to a digraph struct, and two node numbers
 * i and j and returns the value of the change statistic for adding
 * the edge i -- j or arc i -> j (which must not already exist in the graph).
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

/************************* Structural ****************************************/

double changeIsolates(graph_t *g, uint_t i, uint_t j, double lambda);
double changeTwoPath(graph_t *g, uint_t i, uint_t j, double lambda);
double changeLoop(graph_t *g, uint_t i, uint_t j, double lambda);


/************************* Actor attribute (binary) **************************/

double changeInteraction(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);

/********************* Actor attribute (categorical) *************************/

double changeMatching(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeMismatching(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);

/********************* Actor attribute (continuous) *************************/

double changeDiff(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeSum(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);


/********************* Actor attribute (set of categorical) *******************/

double changeJaccardSimilarity(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);



/********************* Dyadic covariate (continuous) *************************/

double changeGeoDistance(graph_t *g, uint_t i, uint_t j);
double changeLogGeoDistance(graph_t *g, uint_t i, uint_t j);
double changeEuclideanDistance(graph_t *g, uint_t i, uint_t j);




/************ Actor attribute interaction (categorical) *********************/

double changeMatchingInteraction(graph_t *g, uint_t i, uint_t j,
                                 uint_t a, uint_t b);



/*************************** Other functions *********************************/

double calcChangeStats(graph_t *g, uint_t i, uint_t j,
                       uint_t n, uint_t n_attr, uint_t n_dyadic,
                       uint_t n_attr_interaction,
                       change_stats_func_t *change_stats_funcs[],
                       double               lambda_values[],
                       attr_change_stats_func_t *attr_change_stats_funcs[],
                       dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                       attr_interaction_change_stats_func_t 
                                        *attr_interaction_change_stats_funcs[],
                       uint_t attr_indices[],
                       double exponent_values[],
                       uint_pair_t attr_interaction_pair_indices[],
                       const double theta[],
                       bool isDelete,
                       double changestats[]);


double jaccard_index(set_elem_e a[], set_elem_e b[], uint_t n);

double *empty_graph_stats(graph_t *g,
                          uint_t n, uint_t n_attr, uint_t n_dyadic,
                          uint_t n_attr_interaction,
                          change_stats_func_t *change_stats_funcs[],
                          double lambda_values[],
                          attr_change_stats_func_t *attr_change_stats_funcs[],
                          dyadic_change_stats_func_t *dyadic_change_stats_funcs[],
                          attr_interaction_change_stats_func_t 
                          *attr_interaction_change_stats_funcs[],
                          uint_t attr_indices[],
                          double exponent_values[],                          
                          uint_pair_t attr_interaction_pair_indices[],
                          double emptystats[]);



#endif /* CHANGESTATISTICSGENERAL_H */

