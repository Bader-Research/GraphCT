#include <stdio.h>
#include <stdlib.h>
#ifdef __MTA__
#include <sys/mta_task.h>
#include <machine/mtaops.h>
#include <machine/runtime.h>
#else
#include "compat/xmt-ops.h"
#endif
#include "defs.h"
#include "globals.h"


int64_t connectedComponents(graph *G)
{
  int64_t count = 0;
  int64_t max_iter = 10;
  int64_t nchanged;

  const int64_t NV = G->numVertices;
  const int64_t NE = G->numEdges;
  const int64_t * restrict sV = G->startVertex;
  const int64_t * restrict eV = G->endVertex;
  int64_t * restrict D = G->marks;

  OMP("omp parallel for")
  for (int64_t i = 0; i < NV; ++i)
    D[i] = i;

  while (max_iter > 0)
  {
    max_iter--;
    nchanged = 0;
    MTA("mta assert nodep")
    OMP("omp parallel for reduction(+:nchanged)")
    for (int64_t k = 0; k < NE; ++k)
    {
      const int64_t i = sV[k];
      const int64_t j = eV[k];
      if (D[i] < D[j])
      {
	D[j] = D[i];
	++nchanged;
      }
    }
    if (!nchanged) break;
  }

  if (nchanged) {
    while (1)
    {
      nchanged = 0;
      MTA("mta assert nodep")
      OMP("omp parallel for reduction(+:nchanged)")
      for (int64_t k = 0; k < NE; ++k)
      {
	const int64_t i = sV[k];
	const int64_t j = eV[k];
	if (D[i] < D[j] && D[j] == D[D[j]])
	{
	  D[D[j]] = D[i];
	  ++nchanged;
	}
      }
      if (!nchanged) break;
      MTA("mta assert nodep")
      OMP("omp parallel for")
      for (int64_t i = 0; i < NV; ++i)
	while (D[i] != D[D[i]])
	  D[i] = D[D[i]];
    }
  }

  MTA("mta assert nodep")
  OMP("omp parallel for reduction(+:count)")
  for (int64_t i = 0; i < NV; ++i)
  {
    while (D[i] != D[D[i]])
      D[i] = D[D[i]];
    if (D[i] == i) ++count;
  }

  return count;
}

