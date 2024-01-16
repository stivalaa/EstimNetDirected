#ifndef CHANGESTATISTICSBIPARTITEUNDIRECTED_H
#define CHANGESTATISTICSBIPARTITEUNDIRECTED_H
/*****************************************************************************
 * 
 * File:    changeStatisticsBipartiteUndirected.h
 * Author:  Alex Stivala
 * Created: May 2022
 *
 * Functions to compute graph change statistics for undirected
 * bipartite (two-mode) graphs. Each function takes a pointer to a
 * graph struct, and two node numbers i and j and returns the value of
 * the change statistic for adding the edge i -- j (which must not
 * already exist in the graph) where i is a node of MODE_A and j is a
 * node of MODE_B.
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
 * Some of these functions are adapted from the original BPNet
 * code by Peng Wang:
 *
 *   Wang P, Robins G, Pattison P. PNet: A program for the simulation
 *   and estimation of exponential random graph models. Melbourne
 *   School of Psychological Science, The University of
 *   Melbourne. 2006. http://www.melnet.org.au/s/PNetManual.pdf
 *
 * And for the definitions of the change statistics:
 *
 *   Wang, P., Sharpe, K., Robins, G. L., & Pattison,
 *   P. E. (2009). Exponential random graph (p∗) models for affiliation
 *   networks. Social Networks, 31(1), 12-25.
 *
 *   Wang, P., Pattison, P., & Robins, G. (2013). Exponential random
 *   graph model specifications for bipartite networks—A dependence
 *   hierarchy. Social networks, 35(2), 211-222.
 * 
 * And also generally:
 * 
 *   Lusher, D., Koskinen, J., & Robins, G. (Eds.). (2013). Exponential
 *   random graph models for social networks: Theory, methods, and
 *   applications. New York, NY: Cambridge University Press.
 * 
 * especially Ch. 10:
 *
 *   Weng, P. (2013). Exponential random graph model extensions:
 *   models for multiple networks and bipartite networks. In
 *   "Exponential random graph models for social networks: Theory,
 *   methods, and applications." (pp. 115-129). New York, NY: Cambridge
 *   University Press.
 *
 * The reference for the MPNet software also defines parameters for
 * bipartite networks as a special case of multilevel networks:
 *
 *   Wang, P., Robins, G., Pattison, P., & Koskinen,
 *   J. (2014). Program for the Simulation and Estimation of (p*)
 *   Exponential Random Graph Models for Multilevel
 *   Networks. Melbourne School of Psychological Sciences, The
 *   University of Melbourne.
 *   https://www.melnet.org.au/s/MPNetManual.pdf
 *
 * As well as the statnet ergm terms, and references for specific
 * change statistics where indicated.
 *
 * Note that in the references above, the two modes are referred to as
 * 'A' and 'P' for 'association' and 'people' respectively, and
 * conventionally represented graphically (red) circle for 'P' nodes
 * and (blue) squares for 'A' nodes. This derives from the common use
 * of bipartite graphs to represent affiliation networks (e.g. company
 * directors on boards) where in the affiliation matrix, the people
 * (P) are rows and their affiliations (A) the columns. However here
 * we are going to use just name the node types A and B, where the
 * type A ndoes correspond to the 'P' in BPNet and type B to the 'A'
 * in BPNet [note this seems to be the other way around in the MPNet
 * naming convention, where e.g. Sa2 in BPNet seems to correspond to
 * XStar2A and Sp2 to XStar2B].
 * So e.g. alternating k-star for for P (k-P-star) is called KSp in BPNet,
 * and changeAltStarsA here (and KSa is changeAltStarsB),
 *
 * b1nodematch and b2nodematch (statnet ergm names) are defined in:
 *
 *    Bomiriya, R. P. (2014). Topics in exponential random graph
 *    modeling. (Doctoral dissertation, Pennsylvania State University).
 *    https://etda.libraries.psu.edu/files/final_submissions/9857
 *
 *    Bomiriya, R. P., Kuvelkar, A. R., Hunter, D. R., & Triebel,
 *    S. (2023). Modeling Homophily in Exponential-Family Random Graph
 *    Models for Bipartite Networks. arXiv preprint
 *    arXiv:2312.05673. https://arxiv.org/abs/2312.05673
 *
 * Do NOT compile with -ffast-math on gcc as we depend on IEEE handling of NaN
 *
 ****************************************************************************/

#include "utils.h"
#include "graph.h"
#include "changeStatisticsTypes.h"

/************************* Utility functions *********************************/
double pow0(uint_t x, double y);

/************************* Structural ****************************************/
double changeBipartiteTwoStarsA(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteTwoStarsB(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteThreeStarsA(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteThreeStarsB(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteAltStarsA(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteAltStarsB(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteAltKCyclesA(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteAltKCyclesB(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteAltK4CyclesA(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteAltK4CyclesB(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteIsolatesA(graph_t *g, uint_t i, uint_t j, double lambda);
double changeBipartiteIsolatesB(graph_t *g, uint_t i, uint_t j, double lambda);

/************************* Actor attribute (binary) **************************/

double changeBipartiteActivityA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteActivityB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);

/*********************** Actor attribute (continuous) ************************/

double changeBipartiteContinuousActivityA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteContinuousActivityB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteTwoPathSumA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteTwoPathSumB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteTwoPathDiffA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteTwoPathDiffB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);


/****************** Actor attribute (continuous, exponent) *******************/

double changeBipartiteDiffBetaA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double alpha);
double changeBipartiteDiffBetaB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double alpha);


/*********************** Actor attribute (categorical) ************************/

double changeBipartiteTwoPathMatchingA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteTwoPathMatchingB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteTwoPathMismatchingA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);
double changeBipartiteTwoPathMismatchingB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double exponent);

/**************** Actor attribute (categorical, exponent) ********************/

double changeBipartiteNodematchAlphaA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double alpha);
double changeBipartiteNodematchAlphaB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double alpha);
double changeBipartiteNodematchBetaA(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double beta);
double changeBipartiteNodematchBetaB(graph_t *g, uint_t i, uint_t j, uint_t a, bool isDelete, double beta);

#endif /* CHANGESTATISTICSBIPARTITEUNDIRECTED_H */
