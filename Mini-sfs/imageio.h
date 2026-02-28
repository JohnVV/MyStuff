/*------------------------------------------------------------------
  imageio.h  –  image I/O utilities
                input  : TIFF  (.tif / .tiff, 8-bit or 16-bit, grayscale
                                or RGB; converted to single-channel float)
                output : TIFF  (.tif, single-channel 32-bit float)

  External dependency: libtiff 4.x  (https://libtiff.gitlab.io/libtiff/)
      Install on Debian/Ubuntu:  sudo apt install libtiff-dev
      Install on Fedora/RHEL  :  sudo dnf install libtiff-devel
  Link with:  -ltiff
------------------------------------------------------------------*/
#pragma once

#include <stdlib.h>
#include <stdio.h>

struct header
{
    int   nrows;
    int   ncols;
    float cellsize;
    float xllcorner;
    float yllcorner;
};

float**         readTIFF (const char* filename, header* hdr);
void            writeTIFF(const char* filename, float** dem, const header* hdr);
float**         fimage   (int nr, int nc);
unsigned char** ucimage  (int nr, int nc);
