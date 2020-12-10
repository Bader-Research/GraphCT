#include "graph_shared.h"
#include "defs.h"
#include "globals.h"
#include "xmt-luc.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netinet/in.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>

void *
shmmap (char * name, int oflags, mode_t mode, int prot, size_t size) 
{
#if !defined(__MTA__)
  int fd = shm_open(name, oflags, mode);
#else
  int fd = open(name, oflags);
#endif

  if(fd == -1) {
    fprintf(stderr, "\nshmmap shm_open error %s\n", strerror(errno)); fflush(stdout);
    return NULL;
  } 

#if !defined(__MTA__)
  int dontcare = ftruncate(fd, size);
  /* silently ignore ftruncate errors */
  
  void * rtn = mmap(NULL, size, prot, MAP_SHARED, fd, 0);
#else
  void * rtn = mmap(NULL, size, prot, MAP_SHARED|MAP_ANON, fd, 0);
#endif
  if(rtn == MAP_FAILED) {
    fprintf(stderr, "\nshmmap mmap error %s\n", strerror(errno)); fflush(stdout);
    return NULL;
  }

  return rtn;
} 

int
shmunmap (char * name, void * ptr, size_t size) 
{
  if(munmap(ptr, size))
    return -1;

#if !defined(__MTA__)
  if(shm_unlink(name))
    return -1;
#else
  if(unlink(name))
    return -1;
#endif

  return 0;
}

void
alloc_shared_graph (graph * G, int64_t NV, int64_t NE, struct graph_shared ** out, char ** loc)
{
  if (!G) return;

  size_t map_size = (3 * NE + 2 * NV + 2) * sizeof (int64_t);

#if !defined(__MTA__)
  sprintf(*loc, "/%lx", (uint64_t)rand());
#else
  char *pwd = xmalloc (sizeof(char) * (MAX_NAME_LEN-16));
  getcwd(pwd, MAX_NAME_LEN-16);
  *loc = xmalloc (sizeof(char) * MAX_NAME_LEN);
  sprintf(*loc, "%s/%lx", pwd, (uint64_t)rand());
#endif
  struct graph_shared * shared = shmmap(*loc, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR,
    PROT_READ | PROT_WRITE, 1 * sizeof(*shared));

  G->numEdges = shared->numEdges = NE;
  G->numVertices = shared->numVertices = NV;
  G->map_size = shared->map_size = map_size;

#if !defined(__MTA__)
  sprintf(shared->startVertex, "/%lx", (uint64_t)rand());
#else
  sprintf(shared->startVertex, "%s/%lx", pwd, (uint64_t)rand());
  free(pwd);
#endif
  G->startVertex = shmmap(shared->startVertex, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR, 
    PROT_READ | PROT_WRITE, map_size);

  G->endVertex = &G->startVertex[NE];
  G->edgeStart = &G->endVertex[NE];
  G->intWeight = &G->edgeStart[NV+2];
  G->marks = xmalloc(NV * sizeof(int64_t));

  *out = shared;
}

void
graph_shared_new (graph * localG, struct graph_shared ** G, char ** out, char * filename, int64_t mode)
{
  int64_t i;
  int64_t NE, NV;
  size_t map_size;
  size_t sz;

  xmt_luc_io_init();
  xmt_luc_stat (filename, &sz);

  if (mode == 4)
  {
    uint32_t * input_graph = (uint32_t *) xmalloc (sz);
    xmt_luc_snapin (filename, input_graph, sz);

    NE = input_graph[0];
    NV = input_graph[1];
    map_size = (3 * NE + 2 * NV + 2) * sizeof (int64_t);

    char *pwd = xmalloc (sizeof(char) * MAX_NAME_LEN);
    getcwd(pwd, MAX_NAME_LEN);
    *out = xmalloc (sizeof(char) * 2 * MAX_NAME_LEN);
    sprintf(*out, "%s/%lx", pwd, (uint64_t)rand());
    struct graph_shared * shared = shmmap(*out, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR,
      PROT_READ | PROT_WRITE, 1 * sizeof(*shared));

    localG->numEdges = shared->numEdges = NE;
    localG->numVertices = shared->numVertices = NV;
    localG->map_size = shared->map_size = map_size;

    sprintf(shared->startVertex, "%s/%lx", pwd, (uint64_t)rand());
    localG->startVertex = shmmap(shared->startVertex, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR, 
      PROT_READ | PROT_WRITE, map_size);

    localG->endVertex = &localG->startVertex[NE];
    localG->edgeStart = &localG->endVertex[NE];
    localG->intWeight = &localG->edgeStart[NV+2];
    localG->marks = xmalloc(NV * sizeof(int64_t));

    uint32_t * start = input_graph + 3;
    uint32_t * eV = input_graph + 3 + NV;
    uint32_t * weight = input_graph + 3 + NV + NE;

    MTA("mta assert nodep")
      for(i = 0; i<NE; i++)
      {
#ifdef __MTA__
	localG->endVertex[i] = (int64_t) ntohl (eV[i]);
	localG->intWeight[i] = (int64_t) ntohl (weight[i]);
#else
	localG->endVertex[i] = (int64_t) eV[i];
	localG->intWeight[i] = (int64_t) weight[i];
#endif
      }

    MTA("mta assert nodep")
      for(i = 0; i<NV; i++)
      {
#ifdef __MTA__
	localG->edgeStart[i] = (int64_t) ntohl (start[i]);
#else
	localG->edgeStart[i] = (int64_t) start[i];
#endif
	localG->marks[i] = 0;
      }

    localG->edgeStart[NV] = NE;

    MTA("mta assert nodep")
      for(i = 0; i<NV; i++) {
	int64_t k;
	int64_t myStart = localG->edgeStart[i];
	int64_t myEnd = localG->edgeStart[i+1];
	MTA("mta assert nodep")
	  for (k = myStart; k < myEnd; ++k)
	    localG->startVertex[k] = i;
      }

    free(pwd);
    free(input_graph);
    *G = shared;
  }
}

void
graph_shared_map (graph * G, struct graph_shared ** shared, char * name)
{
  int64_t NE, NV;
  size_t map_size;

  printf("Connecting to: %s\n", name);
  *shared = shmmap(name, O_RDONLY, S_IRUSR, PROT_READ, 1 * sizeof(struct graph_shared));

  NE = G->numEdges = (*shared)->numEdges;
  NV = G->numVertices = (*shared)->numVertices;
  map_size = G->map_size = (*shared)->map_size;
  G->undirected = (*shared)->undirected;

  G->startVertex = shmmap((*shared)->startVertex, O_RDONLY, S_IRUSR, PROT_READ, map_size);
  
  G->endVertex = &G->startVertex[NE];
  G->edgeStart = &G->endVertex[NE];
  G->intWeight = &G->edgeStart[NV+2];
  G->marks = xmalloc (NV * sizeof(int64_t));
}

void
graph_shared_free (graph *G, struct graph_shared * shared, char * name)
{
  if (!G)
    return;

  int result;
  result = shmunmap (shared->startVertex, G->startVertex, shared->map_size);
  result = shmunmap (name, shared, sizeof(struct graph_shared));
  //free (G->marks);
  //free (G);
  //free (name);
}

void
stinger_shared_unmap (graph *G, struct graph_shared * shared, char * name)
{
}
