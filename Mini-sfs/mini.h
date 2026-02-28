/*------------------------------------------------------------------
  mini.h  –  shared declarations for the shape-from-shading program
------------------------------------------------------------------*/
#ifndef _MINI_H_
#define _MINI_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "imageio.h"

using namespace std;

/* Allow -Ddebug=0 or -Ddebug=1 at compile time; default to 0 */
#ifndef debug
#define debug 0
#endif

/* ----------------------------------------------------------------
   Mathematical constants
---------------------------------------------------------------- */
const double PI  = 3.141592653589793238;
const double R2D = 57.29577952;
const double D2R = 0.017453292;

/* ----------------------------------------------------------------
   Function declarations
---------------------------------------------------------------- */
void usage(void);

void mini_sfs(float** inImg,
              float** initDEM,
              int     nr,
              int     nc,
              double  cellSize,
              double  sunElev,
              double  sunAzim,
              double  lamda,
              int     iterNum);

/* fft2d is now implemented with FFTW 3.3 (fft2d.cpp).
   fftnd has been removed – it was an internal helper of the old
   Numerical Recipes implementation and is no longer needed. */
void fft2d(float* data, int nrows, int ncols, int isign);

#endif  /* _MINI_H_ */
