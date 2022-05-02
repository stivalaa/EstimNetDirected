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
 * the change statistic for adding the edge i -- j where i and j are
 * nodes of differring types (modes), i.e. one is a node of MODE_A and
 * the other of MODE_B.
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
 * Do NOT compile with -ffast-math on gcc as we depend on IEEE handling of NaN
 *
 ****************************************************************************/

#include "utils.h"
#include "graph.h"
#include "changeStatisticsTypes.h"

/************************* Structural ****************************************/

double changeAltStarsA(graph_t *g, uint_t i, uint_t j, double lambda);
double changeAltStarsB(graph_t *g, uint_t i, uint_t j, double lambda);



#endif /* CHANGESTATISTICSBIPARTITEUNDIRECTED_H */

