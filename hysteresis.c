/*------------------------- begin hysteresis.c -------------------------*/
/*******************************************************************************
* FILE: hysteresis.c
* This code was re-written by Mike Heath from original code obtained indirectly
* from Michigan State University. heath@csee.usf.edu (Re-written in 1996).
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0

void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval, int cols);
void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols, float tlow, float thigh, unsigned char *edge);

/*******************************************************************************
* PROCEDURE: follow_edges
* PURPOSE: This procedure traces edges along all paths whose magnitude values
* remain above some specifiable lower threshold.
*******************************************************************************/
void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval, int cols) {
   int i;
   int x[8] = {1,1,0,-1,-1,-1,0,1};
   int y[8] = {0,1,1,1,0,-1,-1,-1};

   for (i = 0; i < 8; i++) {
      unsigned char *tempmapptr = edgemapptr - y[i] * cols + x[i];
      short *tempmagptr = edgemagptr - y[i] * cols + x[i];

      if ((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)) {
         *tempmapptr = (unsigned char) EDGE;
         follow_edges(tempmapptr, tempmagptr, lowval, cols);
      }
   }
}

/*******************************************************************************
* PROCEDURE: apply_hysteresis
* PURPOSE: Finds edges that are above a high threshold or are connected to a high pixel by a path of pixels greater than a low threshold.
*******************************************************************************/
void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols, float tlow, float thigh, unsigned char *edge) {
   int r, c, pos, numedges, highcount, lowthreshold, highthreshold, i, hist[32768];
   short int maximum_mag;

   // Initialize edge map
   for (r = 0, pos = 0; r < rows; r++) {
      for (c = 0; c < cols; c++, pos++) {
         if (nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
         else edge[pos] = NOEDGE;
      }
   }

   // Compute histogram
   for (i = 0; i < 32768; i++) hist[i] = 0;
   for (r = 0, pos = 0; r < rows; r++) {
      for (c = 0; c < cols; c++, pos++) {
         if (edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
      }
   }

   // Find high threshold
   highcount = (int)(numedges * thigh + 0.5);
   for (r = 1, numedges = 0; r < 32768; r++) {
      if (hist[r] != 0) maximum_mag = r;
      numedges += hist[r];
   }
   while ((r < maximum_mag - 1) && (numedges < highcount)) {
      r++;
      numedges += hist[r];
   }
   highthreshold = r;
   lowthreshold = (int)(highthreshold * tlow + 0.5);

   // Trace edges
   for (r = 0, pos = 0; r < rows; r++) {
      for (c = 0; c < cols; c++, pos++) {
         if ((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)) {
            edge[pos] = EDGE;
            follow_edges(edge + pos, mag + pos, lowthreshold, cols);
         }
      }
   }

   // Set all the remaining possible edges to non-edges
   for (r = 0, pos = 0; r < rows; r++) {
      for (c = 0; c < cols; c++, pos++) {
         if (edge[pos] != EDGE) edge[pos] = NOEDGE;
      }
   }
}

/*------------------------- end hysteresis.c -------------------------*/
