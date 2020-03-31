#ifndef IFT_ITERATED_DYNAMIC_OPF_H_
#define IFT_ITERATED_DYNAMIC_OPF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift.h"
#include "iftDataSet.h"

typedef struct ift_graphnode {
    /**  Maximum arc weight from the node to its neighbors */
    float     maxarcw;
    /** Corresponding root node in the graph */
    int       root;
    /** Corresponding training sample in the original dataset */
    int       sample;
    /** List of adjacent nodes */
    iftAdjSet *adj;
    /** List of adjacent nodes on plateaus of density */
    iftSet    *adjplat;
    /** Predecessor node */
    int       pred;     // predecessor node
} iftGraphNode;

typedef struct ift_graph {
    /** Is the graph complete  */
    bool         cpl;
    /** List of nodes in the graph */
    iftGraphNode *node;
    /** List of path value of the nodes */
    float        *pathval;
    /** List of nodes ordered by its path value */
    int          *ordered_nodes;
    /** Number of nodes of the graph */
    int          nnodes;
    /** Priority queue */
    iftFHeap     *Q;
    /** Corresponding dataset */
    iftDataSet   *Z;
} iftGraph;

/**
 * @brief TODO
 * @param Z --- 
 * @param S --- 
 * @param centroids --- 
 * @param k --- 
 * @author David A Cardenas 
 * @date Jan 10th, 2019
 */

void iftUpdateCentroids(iftDataSet* Z, iftDynamicSet **S, int *centroids, int k);

void iftMUpdateCentroids(iftDataSet* Z, iftDynamicSet **S, int *centroids, int k);

/**
 * @brief TODO
 * @param d ---
 * @param f --- 
 * @param n --- 
 * @return TODO
 * @author David A Cardenas 
 * @date Jan 10th, 2019
 */

float iftEuclDistDblFlt(double *d, float *f, int n);

/**
 * @brief TODO
 * @param d ---
 * @param f --- 
 * @param n --- 
 * @return TODO
 * @author David A Cardenas 
 * @date Jan 29th, 2019
 */

float iftEuclDist(float *f1, float *f2, int n);

/**
 * @brief TODO
 * @param S --- 
 * @param Z --- 
 * @param p --- 
 * @author David A Cardenas 
 * @date Jan 10th, 2019
 */

void iftInsertSetDynamicSet(iftDataSet *Z, iftDynamicSet *S, int p);

/**
 * @brief Check the convergence of the clusters' centroids
 * @param centroids --- current values of centroids
 * @param prev_centroids --- previous values of centroids
 * @param k --- number of centroids
 * @return An integer indicating the centroids convergence: (0) no, (1) yes.
 * @author David A Cardenas 
 * @date Jan 10th, 2019
 */

int iftCheckCentroidsConvergence(int *centroids, int *prev_centroids, int k);

/**
 * @brief Cluster a dataset using the Iterated Dynamic OPF algorithm.
 * @param Z --- input dataset
 * @param k --- desired number of clusters
 * @param maxIterations --- maximum number of iterations
 * @return Complete graph.
 * @author David A Cardenas 
 * @date Jan 10th, 2019
 */
  
void *iftIteratedDynamicOPF(iftGraph *graph, int k, int maxIterations);

/**
 * @brief Cluster a dataset using the Iterated Watersheds algorithm.
 * @param Z --- input dataset
 * @param k --- desired number of clusters
 * @param maxIterations --- maximum number of iterations
 * @return Complete graph.
 * @author David A Cardenas 
 * @date Jan 10th, 2019
 */
  
void *iftIteratedWatersheds(iftGraph *graph, int k, int maxIterations);

/**
 * @brief Cluster a dataset using the Iterated Dynamic OPF algorithm (Delaunay).
 * @param graph --- input graph
 * @param k --- desired number of clusters
 * @param maxIterations --- maximum number of iterations
 * @author David A Cardenas 
 * @date Jan 29th, 2019
 */
  
void iftIteratedDynamicOPFDelaunay(iftGraph *graph, int k, int maxIterations);

/**
 * @brief Cluster a dataset using the Iterated Watersheds algorithm (Delaunay).
 * @param graph --- input graph
 * @param k --- desired number of clusters
 * @param maxIterations --- maximum number of iterations
 * @author David A Cardenas 
 * @date Jan 29th, 2019
 */
  
void iftIteratedWatershedsDelaunay(iftGraph *graph, int k, int maxIterations);

/**
 * @brief Create graph.
 * @param Z --- input dataset
 * @return Graph.
 * @author David A Cardenas 
 * @date Jan 30th, 2019
 */

iftGraph *iftCreateGraph(iftDataSet *Z);

/**
 * @brief Destroy graph.
 * @param Z --- input dataset
 * @author David A Cardenas 
 * @date Jan 30th, 2019
 */

void iftDestroyGraph(iftGraph **graph);

/**
 * @brief Set the graph adjacent nodes for Delaunay Triangulation.
 * @param pathname --- delaunay triangulation file
 * @param graph --- input graph
 * @author David A Cardenas 
 * @date Jan 30th, 2019
 */

void iftReadDelaunayTriangulation(const char *pathname, iftGraph *graph);

void iftSetMGraphAdjacencySets(iftGraph *graph, iftMImage *mimg, iftAdjRel *A);

#ifdef __cplusplus
}
#endif

#endif
