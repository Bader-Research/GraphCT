#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "defs.h"
#include "globals.h"
#if !defined(__MTA__)
#include "compat/xmt-ops.h"
#endif
#include "stinger-atomics.h"

int64_t pagerank (graph *G, double * pr, double epsilon, double dampingfactor, int64_t maxiter)
{
  int64_t i;
  int64_t NV = G->numVertices;
  int64_t NE = G->numEdges;
  int64_t * edgeStart = G->edgeStart;
  int64_t * startV = G->startVertex;
  int64_t * endV = G->endVertex;

  int64_t * rev_off = xmalloc ( (NV+2) * sizeof(int64_t));
  int64_t * rev_ind = xmalloc ( NE * sizeof(int64_t));
  int64_t * degree = xmalloc (NV * sizeof(int64_t));
  double * tmp_pr = xmalloc (NV * sizeof(double));

  for (i = 0; i < NV; i++) {
    int64_t myDeg = edgeStart[i+1] - edgeStart[i];
    degree[i] = myDeg;
  }

  for (i = 0; i < NV+2; i++)
    rev_off[i] = 0;

  rev_off+=2;
  for (i = 0; i < NE; i++) {
    stinger_int64_fetch_add(&rev_off[endV[i]], 1);
  }

  for (i = 1; i < NV; i++)
    rev_off[i] += rev_off[i-1];

  rev_off--;

  MTA("mta assert nodep")
  for (i = 0; i < NE; i++) {
    int64_t index = stinger_int64_fetch_add(&rev_off[endV[i]], 1);
    rev_ind[index] = startV[i];
  }

  rev_off--;

  if (rev_off[0] != 0)
    printf("Ahh!\n");

  if (rev_off[NV] != NE)
    printf("Oh no!\n");

  /* Now we have the reversed graph in CSR: rev_off[], ind_off[] */

  for (i = 0; i < NV; i++) {
    pr[i] = 1 / (double) NV;
  }

//  printf("Starting at %lf\n", pr[0]);

  int64_t iter = maxiter;
  double delta = 1;

  while (delta > epsilon && iter > 0) {

    /* calculate my new pagerank based on the ranks of those that point to me */
    for (i = 0; i < NV; i++) {
      int64_t myStart = rev_off[i];
      int64_t myEnd = rev_off[i+1];
      double newPR = 0;

      for (int64_t j = myStart; j < myEnd; j++) {
	int64_t u = rev_ind[j];
	double prU = pr[u];
	int64_t degU = degree[u];
	newPR += (prU / (double) degU);
      }

      tmp_pr[i] = newPR;
    }

    for (i = 0; i < NV; i++) {
      tmp_pr[i] = tmp_pr[i] * dampingfactor + ( (1 - dampingfactor) / (double) NV );
    }

    delta = 0;
    for (i = 0; i < NV; i++) {
      double mydelta = tmp_pr[i] - pr[i];
      if (mydelta < 0)
	mydelta = -mydelta;
      delta += mydelta;
    }

    for (i = 0; i < NV; i++)
      pr[i] = tmp_pr[i];


//    printf("iter: %ld\n", iter);
    iter--;
  }


  free(tmp_pr);
  free(degree);
  free(rev_ind);
  free(rev_off);

  return 0;
}

void SpMV (graph *G, double *vec, double *out, size_t len)
{
  int64_t i;
  int64_t NV = G->numVertices;
  int64_t NE = G->numEdges;
  int64_t * edgeStart = G->edgeStart;
  int64_t * startV = G->startVertex;
  int64_t * endV = G->endVertex;

  for (i = 0; i < NV; i++) {
    out[i] = 0;
  }

  MTA("mta assert nodep")
  for (i = 0; i < NV; i++) {
    int64_t myStart = edgeStart[i];
    int64_t myEnd = edgeStart[i+1];
    double newVal = 0;

    for (int64_t j = myStart; j < myEnd; j++) {
      int64_t myNeighbor = endV[j];
      newVal += vec[myNeighbor];
    }

    out[i] = newVal;
  }

}

int64_t __experimental_pagerank_withSpMV (graph *G, double * pr, double epsilon, double dampingfactor, int64_t maxiter)
{
  int64_t i;
  int64_t NV = G->numVertices;
  int64_t NE = G->numEdges;
  int64_t * edgeStart = G->edgeStart;
  int64_t * startV = G->startVertex;
  int64_t * endV = G->endVertex;

  graph * revGraph;
  alloc_graph(revGraph, NV, NE);
  int64_t * rev_off = revGraph->edgeStart;
  int64_t * rev_ind = revGraph->endVertex;
  int64_t * degree = xmalloc (NV * sizeof(int64_t));
  double * tmp_pr = xmalloc (NV * sizeof(double));

  for (i = 0; i < NV; i++) {
    int64_t myDeg = edgeStart[i+1] - edgeStart[i];
    degree[i] = myDeg;
  }

  for (i = 0; i < NV+2; i++)
    rev_off[i] = 0;

  rev_off+=2;
  for (i = 0; i < NE; i++) {
    stinger_int64_fetch_add(&rev_off[endV[i]], 1);
  }

  for (i = 1; i < NV; i++)
    rev_off[i] += rev_off[i-1];

  rev_off--;

  MTA("mta assert nodep")
  for (i = 0; i < NE; i++) {
    int64_t index = stinger_int64_fetch_add(&rev_off[endV[i]], 1);
    rev_ind[index] = startV[i];
  }

  rev_off--;

  if (rev_off[0] != 0)
    printf("Ahh!\n");

  if (rev_off[NV] != NE)
    printf("Oh no!\n");

  /* Now we have the reversed graph in CSR: rev_off[], ind_off[] */

  for (i = 0; i < NV; i++) {
    pr[i] = 1 / (double) NV;
  }

//  printf("Starting at %lf\n", pr[0]);

  int64_t iter = maxiter;
  double delta = 1;

  while (delta > epsilon && iter > 0) {

    /* calculate my new pagerank based on the ranks of those that point to me */
    SpMV (revGraph, pr, tmp_pr, NV);

    for (i = 0; i < NV; i++) {
      tmp_pr[i] = tmp_pr[i] * dampingfactor + ( (1 - dampingfactor) / (double) NV );
    }

    delta = 0;
    for (i = 0; i < NV; i++) {
      double mydelta = tmp_pr[i] - pr[i];
      if (mydelta < 0)
	mydelta = -mydelta;
      delta += mydelta;
    }

    for (i = 0; i < NV; i++)
      pr[i] = tmp_pr[i];


//    printf("iter: %ld\n", iter);
    iter--;
  }


  free(tmp_pr);
  free(degree);
  free_graph(revGraph); 

  return 0;
}
