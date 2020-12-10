#ifndef  GRAPH_SHARED_H
#define  GRAPH_SHARED_H

#include "defs.h"

#define MAX_NAME_LEN 256

struct graph_shared {
  int64_t numEdges;
  int64_t numVertices;
  char startVertex[2*MAX_NAME_LEN];
  size_t map_size;
  int64_t undirected;
};


void *
shmmap (char * name, int oflags, mode_t mode, int prot, size_t size);

int
shmunmap (char * name, void * ptr, size_t size);

void
graph_shared_new (graph * localG, struct graph_shared ** shared, char ** name, char * filename, int64_t mode);

void
graph_shared_map (graph * G, struct graph_shared ** shared, char * name);

void
graph_shared_free (graph *G, struct graph_shared * shared, char * name);

void
graph_shared_unmap (graph *G, struct graph_shared * shared, char * name);

#endif  /*GRAPH_SHARED_H*/
