#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
// #include <sys/resource.h>
#include <unistd.h>
#include <assert.h>

#ifdef __MTA__
#include <sys/mta_task.h>
#include <machine/runtime.h>
#endif

#include "globals.h"
#include "defs.h"
#include "graph_shared.h"


int main(int argc, char **argv)
{ 
  int64_t i, j, maxI, Vs, NV;

  graph *G;
  struct graph_shared * shared;
  char *name = xmalloc(255 * sizeof(char));

  double maxBC, *BC;

  int64_t *colors;
  int64_t numComp, max, maxV;
  double *trans;
  int64_t diameter = 0;

  /* Temp vars */
  double time, timeBC, timer();
  int64_t scale;

#ifdef __MTA__
  mta_suspend_event_logging();
#endif

  singleKernelOptions options;
  options.cmdsingle = 1;
  options.cmdscript = 0;
  parseCommandLineOptions(argc, argv, &options);

  scale = SCALE = options.scale;


  /*------------------------------------------------------------------------- */
  /*       Preamble -- Untimed                                                */
  /*------------------------------------------------------------------------- */

  /* User Interface: Configurable parameters, and global program control. */
  printf("\nGraphCT - Tools & Kernels for Massive Graph Analysis:\n");
  printf("Running...\n\n");
  fflush (stdout);


  /*------------------------------------------------------------------------- */
  /*     Graph Construction                                                   */
  /*------------------------------------------------------------------------- */

  /* From the input edges, construct the graph 'G'.  */
  printf("\nLoading the graph...\n");
  fflush (stdout);

  time = timer();

  if (options.graph_type == 1)
  {
    G = (graph *) xmalloc(sizeof(graph));
    //graph_shared_new (G, &shared, &name, options.infilename, 4); 
    graphio_b (G, options.infilename, 4, &shared, &name);
    G->undirected = shared->undirected = 1;
  }
  else if (options.graph_type == 5)
  {
    G = (graph *) xmalloc(sizeof(graph));
    graphio_b (G, options.infilename, 8, &shared, &name);
    G->undirected = 1;
  } /*
  else if (options.graph_type == 4)
  {
    G = (graph *) xmalloc(sizeof(graph));
    graphio_el(G, options.infilename);
    G->undirected = 1;
  }
  else if (options.graph_type == 2)
  {
    G = parse_DIMACS_graph(options.infilename);
    G->undirected = 1;
  } */
  else if (options.graph_type == 3)
  {
    char tmp;
    int64_t scale = 0, edgefactor = 0;
    graphSDG SDGdata;
    graph Gdirected;
    sscanf (options.infilename, "%ld%1[.-_: ]%ld", &scale, &tmp, &edgefactor);
    if (scale <= 0 || edgefactor <= 0) {
      fprintf (stderr, "Improper scale:edgefactor : %s\n",
	  options.infilename);
      abort ();
    }
    SCALE = scale;
    getUserParameters (scale, edgefactor);
    genScalData(&SDGdata, 0.55, 0.1, 0.1, 0.25);
    computeGraph (&Gdirected, &SDGdata, &shared, &name);
    free_edgelist (&SDGdata);
    G = makeUndirected (&Gdirected, &shared, &name);
    free_graph (&Gdirected);
  }
  else {
    fprintf (stderr, "Unrecognized graph type.\n");
    abort ();
  }

  graphCheck(G);

  time = timer() - time;

  printf("Time taken is %9.6lf sec.\n", time);
  printf("Graph is available at: %s\n", name);
  printf("\nPress <return> to close the server.\n");
  fflush(stdout);

  getchar();
  printf("\nServer process shutting down...\n");
  graph_shared_free (G, shared, name);
  return 0;
}

