/*------------------------------------------------------------------
  imageio.cpp  –  image I/O utilities
                  input  : TIFF (8- or 16-bit, grayscale or RGB)
                  output : TIFF (single-channel 32-bit float)

  Uses the libtiff 4.x C API (https://libtiff.gitlab.io/libtiff/).
  Build:  g++ ... -ltiff
------------------------------------------------------------------*/

#include "imageio.h"

#include <tiffio.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdexcept>

/* ================================================================
   readTIFF
   Reads an 8-bit or 16-bit TIFF (grayscale or RGB) and returns a
   2-D float** array of luminance values in the range [0, 1].
   Supports:
     • 1-sample  (grayscale) images
     • 3-sample  (RGB) images  – converted via BT.709 coefficients
     • 4-sample  (RGBA) images – alpha channel is ignored
   Contiguous (chunky) and planar (band-separated) configurations
   are both handled.
================================================================ */
float** readTIFF(const char* filename, header* hdr)
{
    TIFF* tif = TIFFOpen(filename, "r");
    if (!tif)
    {
        fprintf(stderr, "readTIFF: cannot open '%s'\n", filename);
        exit(1);
    }

    uint32_t ncols = 0, nrows = 0;
    uint16_t bps   = 8;    /* bits per sample          */
    uint16_t spp   = 1;    /* samples per pixel        */
    uint16_t cfg   = PLANARCONFIG_CONTIG;
    uint16_t photo = PHOTOMETRIC_MINISBLACK;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,      &ncols);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH,     &nrows);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE,   &bps);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetField(tif, TIFFTAG_PLANARCONFIG,    &cfg);
    TIFFGetField(tif, TIFFTAG_PHOTOMETRIC,     &photo);

    if (ncols == 0 || nrows == 0)
    {
        fprintf(stderr, "readTIFF: '%s' has zero dimensions.\n", filename);
        TIFFClose(tif);
        exit(1);
    }
    if (bps != 8 && bps != 16)
    {
        fprintf(stderr,
                "readTIFF: '%s' has unsupported bit depth %u "
                "(only 8 and 16 are supported).\n",
                filename, (unsigned)bps);
        TIFFClose(tif);
        exit(1);
    }
    if (spp > 4)
    {
        fprintf(stderr,
                "readTIFF: '%s' has %u samples/pixel; "
                "only 1–4 are supported.\n",
                filename, (unsigned)spp);
        TIFFClose(tif);
        exit(1);
    }

    hdr->ncols = (int)ncols;
    hdr->nrows = (int)nrows;

    float** img = fimage((int)nrows, (int)ncols);

    /* One scanline buffer, wide enough for spp samples × bytes-per-sample */
    tsize_t scanlineSize = TIFFScanlineSize(tif);
    unsigned char* rowbuf = (unsigned char*)_TIFFmalloc(scanlineSize);
    if (!rowbuf)
    {
        fprintf(stderr, "readTIFF: out of memory for scanline buffer.\n");
        TIFFClose(tif);
        exit(1);
    }

    double maxval = (bps == 16) ? 65535.0 : 255.0;

    for (uint32_t y = 0; y < nrows; y++)
    {
        if (cfg == PLANARCONFIG_CONTIG || spp == 1)
        {
            /* Read all samples for this row in one call */
            if (TIFFReadScanline(tif, rowbuf, y, 0) < 0)
            {
                fprintf(stderr,
                        "readTIFF: error reading row %u from '%s'.\n",
                        (unsigned)y, filename);
                _TIFFfree(rowbuf);
                TIFFClose(tif);
                exit(1);
            }

            for (uint32_t x = 0; x < ncols; x++)
            {
                double r = 0.0, g = 0.0, b = 0.0;

                if (bps == 8)
                {
                    const unsigned char* p = rowbuf + x * spp;
                    if (spp == 1)
                        r = g = b = p[0];
                    else
                    {
                        r = p[0];
                        g = p[1];
                        b = p[2];
                    }
                }
                else /* bps == 16 */
                {
                    const uint16_t* p =
                        reinterpret_cast<const uint16_t*>(rowbuf) + x * spp;
                    if (spp == 1)
                        r = g = b = p[0];
                    else
                    {
                        r = p[0];
                        g = p[1];
                        b = p[2];
                    }
                }

                /* BT.709 luminance */
                img[y][x] = (float)((0.2126 * r + 0.7152 * g + 0.0722 * b)
                                    / maxval);
            }
        }
        else /* PLANARCONFIG_SEPARATE */
        {
            /* Read each plane separately and combine */
            double ch[4] = {0.0, 0.0, 0.0, 0.0};

            for (uint16_t s = 0; s < spp && s < 4; s++)
            {
                if (TIFFReadScanline(tif, rowbuf, y, s) < 0)
                {
                    fprintf(stderr,
                            "readTIFF: error reading row %u plane %u from '%s'.\n",
                            (unsigned)y, (unsigned)s, filename);
                    _TIFFfree(rowbuf);
                    TIFFClose(tif);
                    exit(1);
                }

                for (uint32_t x = 0; x < ncols; x++)
                {
                    double val = (bps == 8)
                        ? (double)rowbuf[x]
                        : (double)reinterpret_cast<uint16_t*>(rowbuf)[x];

                    /* Accumulate into a temporary per-pixel channel array.
                       We need a separate pass, so store the first plane
                       then combine at the end of the plane loop. */
                    (void)val; /* placeholder – see combined loop below */
                }
            }

            /* For planar images, do a proper per-pixel assembly */
            /* Allocate a small scratch buffer for each channel's row */
            unsigned char** planes = (unsigned char**)malloc(spp * sizeof(unsigned char*));
            for (uint16_t s = 0; s < spp; s++)
            {
                planes[s] = (unsigned char*)_TIFFmalloc(scanlineSize);
                if (!planes[s])
                {
                    fprintf(stderr, "readTIFF: out of memory for plane buffer.\n");
                    _TIFFfree(rowbuf);
                    TIFFClose(tif);
                    exit(1);
                }
                TIFFReadScanline(tif, planes[s], y, s);
            }

            for (uint32_t x = 0; x < ncols; x++)
            {
                double r = 0.0, g = 0.0, b = 0.0;
                if (bps == 8)
                {
                    r = planes[0][x];
                    g = (spp > 1) ? planes[1][x] : r;
                    b = (spp > 2) ? planes[2][x] : r;
                }
                else
                {
                    r = reinterpret_cast<uint16_t*>(planes[0])[x];
                    g = (spp > 1) ? reinterpret_cast<uint16_t*>(planes[1])[x] : r;
                    b = (spp > 2) ? reinterpret_cast<uint16_t*>(planes[2])[x] : r;
                }
                img[y][x] = (float)((0.2126 * r + 0.7152 * g + 0.0722 * b)
                                    / maxval);
            }

            for (uint16_t s = 0; s < spp; s++)
                _TIFFfree(planes[s]);
            free(planes);
        }
    }

    _TIFFfree(rowbuf);
    TIFFClose(tif);

    printf("Read TIFF: %s  (%u x %u, %u bps, %u spp)\n",
           filename, ncols, nrows, (unsigned)bps, (unsigned)spp);
    return img;
}

/* ================================================================
   writeTIFF
   Writes a single-channel 32-bit IEEE float TIFF.
   The channel is tagged PHOTOMETRIC_MINISBLACK with SAMPLEFORMAT_IEEEFP.
================================================================ */
void writeTIFF(const char* filename, float** dem, const header* hdr)
{
    int nrows = hdr->nrows;
    int ncols = hdr->ncols;

    TIFF* tif = TIFFOpen(filename, "w");
    if (!tif)
    {
        fprintf(stderr, "writeTIFF: cannot open '%s' for writing.\n", filename);
        exit(1);
    }

    /* Basic geometry */
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,      (uint32_t)ncols);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH,     (uint32_t)nrows);

    /* Single-channel 32-bit float */
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, (uint16_t)1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,   (uint16_t)32);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,    SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,     PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG,    PLANARCONFIG_CONTIG);

    /* Row-by-row (no compression for maximum compatibility) */
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP,    (uint32_t)1);
    TIFFSetField(tif, TIFFTAG_COMPRESSION,     COMPRESSION_NONE);

    /* Optional spatial metadata from the header */
    if (hdr->cellsize > 0.0f)
    {
        TIFFSetField(tif, TIFFTAG_XRESOLUTION,  (float)(1.0f / hdr->cellsize));
        TIFFSetField(tif, TIFFTAG_YRESOLUTION,  (float)(1.0f / hdr->cellsize));
        TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);
    }

    for (int y = 0; y < nrows; y++)
    {
        if (TIFFWriteScanline(tif, dem[y], (uint32_t)y, 0) < 0)
        {
            fprintf(stderr,
                    "writeTIFF: error writing row %d to '%s'.\n", y, filename);
            TIFFClose(tif);
            exit(1);
        }
    }

    TIFFClose(tif);

    printf("Wrote TIFF: %s  (%d x %d, float32)\n", filename, ncols, nrows);
}

/* ================================================================
   Memory allocation helpers  (unchanged)
================================================================ */

float** fimage(int nr, int nc)
{
    float** x = (float**)malloc((size_t)nr * sizeof(float*));
    if (!x)
    {
        fprintf(stderr, "fimage: out of memory (row pointers).\n");
        exit(1);
    }
    x[0] = (float*)malloc((size_t)nr * nc * sizeof(float));
    if (!x[0])
    {
        fprintf(stderr, "fimage: out of memory (data block).\n");
        exit(1);
    }
    for (int i = 1; i < nr; i++)
        x[i] = x[i - 1] + nc;
    return x;
}

unsigned char** ucimage(int nr, int nc)
{
    unsigned char** x =
        (unsigned char**)malloc((size_t)nr * sizeof(unsigned char*));
    if (!x)
    {
        fprintf(stderr, "ucimage: out of memory (row pointers).\n");
        exit(1);
    }
    x[0] = (unsigned char*)malloc((size_t)nr * nc * sizeof(unsigned char));
    if (!x[0])
    {
        fprintf(stderr, "ucimage: out of memory (data block).\n");
        exit(1);
    }
    for (int i = 1; i < nr; i++)
        x[i] = x[i - 1] + nc;
    return x;
}
