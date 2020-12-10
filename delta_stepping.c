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

#define DEBUG 1
#define INIT_SIZE 100
#define FREEMEM 1

static double INFTY;

static int64_t
relax (int64_t * R,
       double * dR,
       int64_t * RLoc,
       int64_t Rcount,
       graph * G,
       double * d,
       int64_t ** B,
       int64_t currBucketNum,
       int64_t numBuckets,
       int64_t * Bsize,
       int64_t * Bcount,
       int64_t * Braise,
       int64_t * Bdel,
       int64_t * vBucketMap,
       int64_t * vBucketLoc,
       double delta);

int64_t delta_stepping (graph * G, double * d, int64_t s, double * cs)
{
  double delta;
  int64_t i, j, m, n;
  int64_t *numEdges, *endV, *intWeight;
  double *W, *dR;
  int64_t **B, numBuckets, *Bcount, *Bsize, *Braise, *Bdel, *Bt, *Btemp,
      Btemp_count, Bdel_t, Bcount_t;
  int64_t currBucketNum, lastBucketNum, *vBucketMap, *vBucketLoc;
  int64_t *S, *SLoc, Scount, *R, *RLoc, Rcount;
  int64_t *memBlock;
  double delta_l, INFTY_l;
  FILE* outfp;

  int64_t not_connected;

#if DEBUG
  int64_t numPhases, *numVertsInPhase, numRelaxations, numRelaxationsHeavy, maxBucketNum;
#endif

  double Wsum, checksum;

  m = G->numEdges;
  n = G->numVertices;

  numEdges = G->edgeStart;
  endV  = G->endVertex;
  intWeight = G->intWeight;
  W = xmalloc (m * sizeof(*W));

  /* Convert integer edge weights in the graph to doubles */
  int64_t maxWt = 0;
  for (i = 0; i < m; i++)
  {
    if (intWeight[i] > maxWt)
      maxWt = intWeight[i];
  }

  MTA("mta assert nodep")
  for (i = 0; i < m; i++)
  {
    W[i] = ((double) intWeight[i]/maxWt);
  }


  INFTY = ((double) m) + 1.0;
  // delta = ((double) n)/((double) 2*m);
  // delta = 0.3;
  delta = *cs;
  delta_l = delta;
  INFTY_l = INFTY;
  // numBuckets = 30;
  // numBuckets = 50000;
  numBuckets = (int64_t) ((double) (1.0/delta) * ((double)n/2));
  // numBuckets = (int64_t) ((1.0/delta)) * 3 * ((int64_t) ceil(log10(n)/log10(2.0)));

#if DEBUG
  fprintf(stderr, "source: %ld\n", s);
  fprintf(stderr, "delta: %lf\n", delta);
  fprintf(stderr, "No. of buckets: %ld\n", numBuckets);

  /* Assuming a max of n phases */
  /* Storing the no. of verts visited in each phase 
     for statistics */
  numVertsInPhase = xcalloc (n, sizeof(int64_t));
#endif


  /* Memory allocation */

  B = xmalloc (numBuckets*sizeof(int64_t *));
  Bsize = xcalloc (numBuckets, sizeof(int64_t));
  Bcount = xcalloc (numBuckets, sizeof(int64_t));
  Braise = xcalloc (4*numBuckets, sizeof(int64_t));
  Bdel = xcalloc (numBuckets, sizeof(int64_t));

#ifdef __MTA__
  memBlock = xmalloc ((9*n+1)*sizeof(int64_t));
#else
  if (sizeof(int64_t) == 4) 
    memBlock = xmalloc ((11*n+2)*sizeof(int64_t));
  if (sizeof(int64_t) == 8)
    memBlock = xmalloc ((9*n+1)*sizeof(int64_t));
#endif

  S = memBlock;
  R = memBlock + n;
  SLoc = memBlock + 2*n;
  RLoc = memBlock + 3*n;
  Btemp = memBlock + 4*n;
  vBucketMap = memBlock + 5*n;
  vBucketLoc = memBlock + 6*n;
  dR = (double *) (memBlock + 7*n);


  /* Initializing the bucket data structure */
  for (i = 0; i < n; i++)
  {
    vBucketMap[i] = -1;
    d[i] = INFTY_l;   
    RLoc[i] = -1;
    SLoc[i] = -1;
  }

  d[n] = INFTY_l;


  R[0] = s;
  dR[0] = 0;
  Rcount = 0;
  Scount = 0; 

  lastBucketNum = relax(R, dR, RLoc, 1, G, d, B, 0, numBuckets, Bsize, Bcount, Braise, Bdel, vBucketMap, vBucketLoc, delta);

#if DEBUG
  numRelaxations = 1;
  numRelaxationsHeavy = 0;    
  numPhases = 0;
  maxBucketNum = 0;
#endif

  currBucketNum = 0;
  while (currBucketNum <= lastBucketNum) {

    if (Bcount[currBucketNum] == 0) {
      currBucketNum++;
      continue;
    }

    /* Handle light edges */
    while (Bcount[currBucketNum] != 0) {

      Bcount_t = Bcount[currBucketNum];
      Bdel_t = Bdel[currBucketNum];
      Bt = B[currBucketNum];
      Btemp_count = 0;
      Rcount = 0;

      if (Bdel_t == Bcount_t) {
	Btemp_count = 0;

	/* The bucket doesn't have a lot of empty spots */
      } else if (Bdel_t < Bcount_t/3 + 2) {

	Btemp_count = Bcount_t;

	if (Bcount_t > 30)
	{
	  MTA("mta assert parallel")
	  for (i = 0; i < Bcount_t; i++) {
	    double du;
	    int64_t u, v, rlv, start, pos, end;
	    u = Bt[i];
	    du = d[u];
	    start = numEdges[u];
	    end = numEdges[u+1];

	    MTA("mta assert nodep")
	    for (j = start; j < end; j++) {
	      v = endV[j];
	      if (du + W[j] < d[v]) {
		rlv = readfe(&RLoc[v]);
		if (rlv == -1) {
		  pos = stinger_int64_fetch_add(&Rcount, 1);
		  R[pos] = v;
		  dR[pos] = du + W[j];
		  writeef(&RLoc[v], pos);
		} else {
		  if (du + W[j] < dR[rlv])
		    dR[rlv] = du + W[j];
		  writeef(&RLoc[v], rlv);
		}
	      }
	    }
	  }

	} else {
	  for (i = 0; i < Bcount_t; i++) {
	    double du;
	    int64_t u, v, rlv, start, pos, end;
	    u = Bt[i];
	    du = d[u];
	    start = numEdges[u];
	    end = numEdges[u+1];

	    MTA("mta assert nodep")
	    for (j = start; j < end; j++) {
	      v = endV[j];
	      if (du + W[j] < d[v]) {
		rlv = readfe(&RLoc[v]);
		if (rlv == -1) {
		  pos = stinger_int64_fetch_add(&Rcount, 1);
		  R[pos] = v;
		  dR[pos] = du + W[j];
		  writeef(&RLoc[v], pos);
		} else {
		  if (du + W[j] < dR[rlv])
		    dR[rlv] = du + W[j];
		  writeef(&RLoc[v], rlv);
		}
	      }
	    }
	  }
	}

	/* Add to S */
	if (Bcount_t > 30)
	{
	  MTA("mta assert parallel")
	  for (i = 0; i < Bcount_t; i++) {
	    int64_t slv;
	    int64_t Gn = n;
	    int64_t u = Bt[i];

	    if (u == Gn) 
	      continue;

	    slv = readfe(&SLoc[u]);

	    /* Add the vertex to S */
	    if (slv == -1) {
	      int64_t pos = stinger_int_fetch_add(&Scount, 1);
	      S[pos] = u;
	      writeef(&SLoc[u], pos);
	    } else {
	      writeef(&SLoc[u], slv);
	    }

	    vBucketMap[u] = -1;
	    vBucketLoc[u] = -1;
	  }

	} else {

	  for (i = 0; i < Bcount_t; i++) {
	    int64_t slv;
	    int64_t Gn = n;
	    int64_t u = Bt[i];

	    if (u == Gn) 
	      continue;

	    slv = readfe(&SLoc[u]);

	    /* Add the vertex to S */
	    if (slv == -1) {
	      int64_t pos = stinger_int64_fetch_add(&Scount, 1);
	      S[pos] = u;
	      writeef(&SLoc[u], pos);
	    } else {
	      writeef(&SLoc[u], slv);
	    }

	    vBucketMap[u] = -1;
	    vBucketLoc[u] = -1;
	  }
	}

      } else {        

	/* Bdel_t > Bcount_t/3  */
	/* There are a significant number of empty spots in the bucket.
	 * So we get the non-empty vertices and store them in a compact
	 * array */

	int64_t Gn = n;

	if (Bcount_t > 30) {

	  MTA("mta assert nodep")
	  MTA("mta interleave schedule")
	  for (i = 0; i < Bcount_t; i++) {
	    int64_t u = Bt[i];
	    if (u != Gn) {
	      int64_t pos = stinger_int64_fetch_add(&Btemp_count, 1);
	      Btemp[pos] = u;
	    }
	  }

	} else {

	  for (i = 0; i < Bcount_t; i++) {
	    int64_t u = Bt[i];
	    if (u != Gn) {
	      int64_t pos = stinger_int64_fetch_add(&Btemp_count, 1);
	      Btemp[pos] = u;
	    }
	  }
	} 

	/* The nested loop can be collapsed, but this results
	 * in a lot of hotspots */

	if (Btemp_count > 30) {

	  MTA("mta assert parallel")
	  for (i = 0; i < Btemp_count; i++) {
	    int64_t u = Btemp[i];
	    double du = d[u];
	    int64_t start, end;
	    start = numEdges[u];
	    end = numEdges[u+1];

	    for (j = start; j < end; j++) {
	      int64_t v = endV[j];
	      if (du + W[j] < d[v]) {
		int64_t rlv = readfe(&RLoc[v]);
		if (rlv == -1) {
		  int64_t pos = stinger_int64_fetch_add(&Rcount, 1);
		  R[pos] = v;
		  dR[pos] = du + W[j];
		  writeef(&RLoc[v], pos);
		} else {
		  if (du + W[j] < dR[rlv])
		    dR[rlv] = du + W[j];
		  writeef(&RLoc[v], rlv);
		}
	      }
	    }   
	  }

	} else {

	  for (i = 0; i < Btemp_count; i++) {

	    int64_t u = Btemp[i];
	    double du = d[u];
	    int64_t start, end;

	    start = numEdges[u];
	    end = numEdges[u+1];

	    for (j = start; j < end; j++) {
	      int64_t v = endV[j];
	      if (du + W[j] < d[v]) {
		int64_t rlv = readfe(&RLoc[v]);
		if (rlv == -1) {
		  int64_t pos = stinger_int64_fetch_add(&Rcount, 1);
		  R[pos] = v;
		  dR[pos] = du + W[j];
		  writeef(&RLoc[v], pos);
		} else {
		  if (du + W[j] < dR[rlv])
		    dR[rlv] = du + W[j];
		  writeef(&RLoc[v], rlv);
		}
	      }
	    }   
	  }
	}

	if (Btemp_count > 30)
	{     
	  MTA("mta assert parallel")
	  for (i = 0; i < Btemp_count; i++) {

	    int64_t slv;
	    int64_t u = Btemp[i];
	    slv = readfe(&SLoc[u]);

	    /* Add the vertex to S */   

	    if (slv == -1) {
	      int64_t pos = stinger_int64_fetch_add(&Scount, 1);
	      S[pos] = u;
	      writeef(&SLoc[u], pos); 
	    } else {
	      writeef(&SLoc[u], slv);
	    }

	    vBucketMap[u] = -1;
	    vBucketLoc[u] = -1;
	  }

	} else {

	  for (i = 0; i < Btemp_count; i++) {

	    int64_t slv;
	    int64_t u = Btemp[i];
	    slv = readfe(&SLoc[u]);

	    /* Add the vertex to S */   

	    if (slv == -1) {
	      int64_t pos = stinger_int64_fetch_add(&Scount, 1);
	      S[pos] = u;
	      writeef(&SLoc[u], pos); 
	    } else {
	      writeef(&SLoc[u], slv);
	    }

	    vBucketMap[u] = -1;
	    vBucketLoc[u] = -1;
	  }
	}
      } /* end of if .. then .. else */
      /* We have collected all the light edges in R */

#if DEBUG
      if (Btemp_count != 0)
	maxBucketNum = currBucketNum;
#endif

      Bcount[currBucketNum] = 0;
      Bdel[currBucketNum] = 0;

#if DEBUG
      numVertsInPhase[numPhases++] = Rcount;
      numRelaxations += Rcount;
#endif
      /* Relax all light edges */
      if (Rcount != 0) 
	lastBucketNum = relax(R, dR, RLoc, Rcount, G, d, B, currBucketNum, numBuckets, Bsize, Bcount, Braise, Bdel, vBucketMap, vBucketLoc, delta);
    }       

    Rcount = 0;
    /* Collect heavy edges int64_to R */
    if (Scount > 10)
    {
      MTA("mta assert parallel")
      for (i = 0; i < Scount; i++) {
	int64_t u = S[i];
	double du = d[u];   
	int64_t start = numEdges[u];
	int64_t end = numEdges[u+1];
	SLoc[u] = -1;

	for (j = start; j < end; j++) {

	  int64_t v = endV[j];
	  if (W[j] + du < d[v]) {
	    int64_t rlv = readfe(&RLoc[v]);
	    if (rlv == -1) {
	      int64_t pos = stinger_int64_fetch_add(&Rcount, 1);
	      R[pos] = endV[j];
	      dR[pos] = d[u] + W[j];
	      writeef(&RLoc[v], pos);
	    } else {
	      if (du + W[j] < dR[rlv])
		dR[rlv] = du + W[j];
	      writeef(&RLoc[v], rlv);
	    }
	  }
	}
      }

    } else {

      for (i = 0; i < Scount; i++) {
	int64_t u = S[i];
	double du = d[u];   
	int64_t start = numEdges[u];
	int64_t end = numEdges[u+1];
	SLoc[u] = -1;

	for (j = start; j < end; j++) {

	  int64_t v = endV[j];
	  if (W[j] + du < d[v]) {
	    int64_t rlv = readfe(&RLoc[v]);
	    if (rlv == -1) {
	      int64_t pos = stinger_int64_fetch_add(&Rcount, 1);
	      R[pos] = endV[j];
	      dR[pos] = d[u] + W[j];
	      writeef(&RLoc[v], pos);
	    } else {
	      if (du + W[j] < dR[rlv])
		dR[rlv] = du + W[j];
	      writeef(&RLoc[v], rlv);
	    }
	  }
	}
      }
    }

    Scount = 0; 
    /* Relax heavy edges */
    if (Rcount != 0) {
      lastBucketNum = relax(R, dR, RLoc, Rcount, G, d, B, currBucketNum, numBuckets, Bsize, Bcount, Braise, Bdel, vBucketMap, vBucketLoc, delta);
    } 

#if DEBUG
    numRelaxationsHeavy += Rcount;
#endif
  }


  /* Compute d checksum */
  checksum = 0;
  not_connected = 0; 

  MTA("mta assert parallel")
  for (i = 0; i < n; i++) {
    if (d[i] < INFTY_l)
      checksum = checksum + d[i];
    else
      stinger_int64_fetch_add(&not_connected, 1);
  }
  *cs = checksum;

#if DEBUG
  /* Compute W checksum */
  Wsum = 0;
  for (i = 0; i < m; i++) {
    Wsum = Wsum + W[i];
  }

  *cs = checksum;
  fprintf(stderr, "d checksum: %lf, W checksum: %lf, Avg. distance %lf\n",
      checksum, Wsum, checksum/(n-not_connected));

  fprintf(stderr, "Last non-empty bucket: %ld\n", maxBucketNum);
  fprintf(stderr, "No. of phases: %ld\n", numPhases);

  fprintf(stderr, "No. of light relaxations: %ld\n", numRelaxations);
  fprintf(stderr, "No. of heavy relaxations: %ld\n", numRelaxationsHeavy); 
  fprintf(stderr, "Avg. no. of light edges relaxed in a phase: %lf\n",
      (double) numRelaxations/numPhases);
  fprintf(stderr, "Avg. no. of heavy edges relaxed per bucket: %lf\n",
      (double) numRelaxationsHeavy/(maxBucketNum+1));
  fprintf(stderr, "Total no. of relaxations: %ld\n\n", numRelaxations+numRelaxationsHeavy);
#endif 

#if FREEMEM 
  /* Free memory */    
  MTA("mta assert parallel")
  for (i = 0; i < numBuckets; i++) {
    if (Bsize[i] != 0)
      free(B[i]);
  }

  free(W);
  free(B);
  free(Bsize);
  free(Bcount);
  free(Braise); 
  free(Bdel);
  free(memBlock);

/*
#if DEBUG
  outfp = fopen("phase_plot.txt", "w");
  for (i = 0; i < numPhases; i++) {
    if (numVertsInPhase[i] > 10)
      fprintf(outfp, "%ld\t%ld\n", i, numVertsInPhase[i]);
  }
  fclose(outfp);
  free(numVertsInPhase);
#endif
*/
#endif

  return not_connected;    
}


MTA("mta inline")
static int64_t
relax (int64_t * R,
       double * dR,
       int64_t * RLoc,
       int64_t Rcount,
       graph * G,
       double * d,
       int64_t ** B,
       int64_t currBucketNum,
       int64_t numBuckets,
       int64_t * Bsize,
       int64_t * Bcount,
       int64_t * Braise,
       int64_t * Bdel,
       int64_t * vBucketMap,
       int64_t * vBucketLoc,
       double delta)
{
  int64_t i;
  double delta_l = delta;
  double INFTY_l = INFTY;
  int64_t Gn = G->numVertices;
  int64_t lastBucketNum = -1;

  if (Rcount > 30)
  {
    MTA("mta assert nodep")
    MTA("mta interleave schedule")
    for (i = 0; i < Rcount; i++) {
      int64_t v = R[i];
      int64_t bn = (int64_t) floor(dR[i]/delta_l);
      int64_t bn_old = vBucketMap[v];
      int64_t offset = numBuckets * (i & 3);

      if (bn > lastBucketNum)
	lastBucketNum = bn;
      /*
	 if (bn >= numBuckets) {
	 fprintf(stderr, "Error: relaxation failed, bn: %d, numBuckets: %d\n", bn, numBuckets); 
	 exit(1); 
	 }
       */
      RLoc[v] = (i & 3) * Gn + stinger_int64_fetch_add(&Braise[bn+offset], 1);

      if ((d[v] < INFTY_l) && (bn_old != -1)) {
	B[bn_old][vBucketLoc[v]] = Gn;
	Bdel[bn_old]++;
      }
    }
  } else {
    for (i = 0; i < Rcount; i++) {
      int64_t v = R[i];
      int64_t bn = (int64_t) floor(dR[i]/delta_l);
      int64_t bn_old = vBucketMap[v];
      int64_t offset = numBuckets * (i & 3);

      if (bn > lastBucketNum)
	lastBucketNum = bn;
      /*
	 if (bn >= numBuckets) {
	 fprintf(stderr, "Error: relaxation failed, bn: %d, numBuckets: %d\n", bn, numBuckets); 
	 exit(1); 
	 }
       */
      RLoc[v] = (i & 3) * Gn + stinger_int64_fetch_add(&Braise[bn+offset], 1);

      if ((d[v] < INFTY_l) && (bn_old != -1)) {
	B[bn_old][vBucketLoc[v]] = Gn;
	Bdel[bn_old]++;
      }
    }
  }

  lastBucketNum++;
  // fprintf(stderr, "[%d %d] ", currBucketNum, lastBucketNum);
  for (i = currBucketNum; i < lastBucketNum; i++) {
    int64_t *Bi = B[i];
    int64_t Bsize_i = Bsize[i];
    int64_t size_incr = Braise[i] + Braise[i+numBuckets] + 
      Braise[i+2*numBuckets] + Braise[i+3*numBuckets];

    if ((size_incr > 0) && (Bcount[i] + size_incr >= Bsize[i])) {

      int64_t Bsize_tmp = Bcount[i] + size_incr + INIT_SIZE;
      int64_t* Bt = xmalloc (Bsize_tmp*sizeof(int64_t));

      if (Bsize_i != 0) {
	int64_t j;
	if (Bsize_i > 30) {
	  MTA("mta assert nodep")
	  MTA("mta interleave schedule")
	  for (j = 0; j < Bsize_i; j++)
	    Bt[j] = Bi[j];
	} else {
	  for (j = 0; j < Bsize_i; j++)
	    Bt[j] = Bi[j];
	}
	free(Bi);
      }

      B[i] = Bt;
      Bsize[i] = Bsize_tmp;
    }
  }

  if (Rcount > 30)
  { 
    MTA("mta assert nodep")
    MTA("mta interleave schedule")
    for (i = 0; i < Rcount; i++) {
      int64_t v = R[i];
      double x = dR[i];
      int64_t loc = RLoc[v];
      int64_t locDiv = loc/Gn;
      int64_t locMod = loc % Gn;
      int64_t bn = (int64_t) floor(x/delta_l);
      int64_t pos = Bcount[bn] + locMod;
      pos += (locDiv >= 1) * Braise[bn];
      pos += (locDiv >= 2) * Braise[bn+numBuckets];
      pos += (locDiv >= 3) * Braise[bn+2*numBuckets]; 

      B[bn][pos] = v;
      vBucketLoc[v] = pos;
      vBucketMap[v] = bn;
      d[v] = x;
      RLoc[v] = -1;
    }
  } else {
    for (i = 0; i < Rcount; i++) {
      int64_t v = R[i];
      double x = dR[i];
      int64_t loc = RLoc[v];
      int64_t locDiv = loc/Gn;
      int64_t locMod = loc % Gn;
      int64_t bn = (int64_t) floor(x/delta_l);
      int64_t pos = Bcount[bn] + locMod;
      pos += (locDiv >= 1) * Braise[bn];
      pos += (locDiv >= 2) * Braise[bn+numBuckets];
      pos += (locDiv >= 3) * Braise[bn+2*numBuckets]; 

      B[bn][pos] = v;
      vBucketLoc[v] = pos;
      vBucketMap[v] = bn;
      d[v] = x;
      RLoc[v] = -1;
    }
  }

  if (lastBucketNum - currBucketNum > 30)
  {
    MTA("mta parallel single processor")
    MTA("mta assert nodep")
    MTA("mta interleave schedule")
    for (i = currBucketNum; i < lastBucketNum; i++) {
      Bcount[i] += Braise[i] + Braise[i+numBuckets] + Braise[i+2*numBuckets] +
	Braise[i+3*numBuckets];

      Braise[i] = 0;
      Braise[i+numBuckets] = 0;
      Braise[i+2*numBuckets] = 0;
      Braise[i+3*numBuckets] = 0; 
    }
  } else {

    for (i = currBucketNum; i < lastBucketNum; i++) {
      Bcount[i] += Braise[i] + Braise[i+numBuckets] + Braise[i+2*numBuckets] +
	Braise[i+3*numBuckets];

      Braise[i] = 0;
      Braise[i+numBuckets] = 0;
      Braise[i+2*numBuckets] = 0;
      Braise[i+3*numBuckets] = 0; 
    }
  }

  return lastBucketNum;

}
