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

int64_t bellman_ford (graph *G, int64_t * d, int64_t s, int64_t * cs)
{
  int64_t i;
  int64_t NV = G->numVertices;
  int64_t NE = G->numEdges;
  int64_t * startV = G->startVertex;
  int64_t * endV = G->endVertex;
  int64_t * W = G->intWeight;
  int64_t INFTY;

  fprintf(stderr, "source: %ld\n", s);

  INFTY = INT64_MAX;
  for (i = 0; i < NV; i++) {
    d[i] = INFTY;
  } 
  d[s] = 0;

  for (i = 0; i < NV-1; i++) {
    int64_t j; 

    MTA("mta assert nodep")
    MTA("mta interleave schedule")
    for (j = 0; j < NE; j++) {
      int64_t u = startV[j];
      int64_t v = endV[j];
      int64_t du = d[u];
      int64_t w = W[j];
      int64_t dv = readfe(&d[v]);
      if (du + w < dv) 
	writeef(&d[v], du + w);
      else
	writeef(&d[v], dv);
    } 
  }

  /* Compute d checksum */
  int64_t checksum = 0;
  int64_t not_connected = 0; 

  MTA("mta assert parallel")
  for (i = 0; i < NV; i++) {
    if (d[i] < INFTY)
      checksum = checksum + d[i];
    else
      stinger_int64_fetch_add(&not_connected, 1);
  }
  *cs = checksum;

  return not_connected;    
}
