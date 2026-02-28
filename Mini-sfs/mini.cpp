/****************************************************************************
NAME:    Mini_sfs  –  Shape-from-shading via iterative minimisation

SYNOPSIS:
  mini [options] <input.tif> <output.tif>

  Required arguments (positional):
    input.tif          Input image (TIFF, 8- or 16-bit, grayscale or RGB)
    output.tif         Output elevation model (TIFF, single-channel float32)

  Required named arguments:
    -e <sunElevAngle>  Sun elevation angle in decimal degrees
    -a <sunAzimAngle>  Sun azimuth angle  in decimal degrees (N=0, clockwise)
    -l <lambda>        Regularisation parameter (typical range 0.4 – 50)
    -n <iterNum>       Maximum number of iterations (e.g. 50)

  Optional:
    -d <initDem.tif>   Initial DEM as a TIFF image (grayscale/luminance)
    -c <cellsize>      Ground-sample distance in metres per pixel (default 1.0)
    -x <xllcorner>     Easting  of lower-left corner (default 0.0)
    -y <yllcorner>     Northing of lower-left corner (default 0.0)

EXAMPLES:
  Mini_sfs -e 19.3 -a 287.2 -l 1.5 -n 50 spot128.tif ts_dem.tif
  Mini_sfs -e 19.3 -a 287.2 -l 1.5 -n 50 -d init.tif -c 30.0 spot128.tif ts_dem.tif

ALGORITHM REFERENCE:
  Liu, H. "Derivation of surface topography and terrain parameters from a
  single satellite image using shape-from-shading" – Computers & Geosciences : year-2002

BUILD:
  g++ -O2 -Ddebug=0 \
      mini.cpp mini_sfs.cpp fft2d.cpp imageio.cpp \
      -I. -lm -ltiff -lfftw3f_threads -lfftw3f -lpthread -o Mini_sfs
****************************************************************************/

#include "mini.h"
#include <string.h>

/* ----------------------------------------------------------------- */
void usage(void)
{
    cerr << "\nUsage:\n"
         << "  mini [options] <input.tif> <output.tif>\n\n"
         << "Required options:\n"
         << "  -e <deg>   Sun elevation angle (decimal degrees)\n"
         << "  -a <deg>   Sun azimuth  angle  (decimal degrees, N=0 CW)\n"
         << "  -l <val>   Lambda – regularisation parameter (0.4 – 50)\n"
         << "  -n <int>   Maximum iteration count (e.g. 50)\n\n"
         << "Optional:\n"
         << "  -d <file>  Initial DEM as TIFF (default: flat zero surface)\n"
         << "  -c <val>   Cell size in metres per pixel (default 1.0)\n"
         << "  -x <val>   Easting  of lower-left corner (default 0.0)\n"
         << "  -y <val>   Northing of lower-left corner (default 0.0)\n\n"
         << "Examples:\n"
         << "  mini -e 19.3 -a 287.2 -l 1.5 -n 50 spot128.tif dem.tif\n"
         << "  mini -e 45.0 -a 135.0 -l 2.0 -n 100 -d init.tif -c 30 "
            "image.tif result.tif\n\n";
    exit(1);
}

/* ================================================================= */
int main(int argc, char** argv)
{
    /* ---- default values --------------------------------------- */
    char   inName[256]      = {0};
    char   outName[256]     = {0};
    char   initDemName[256] = {0};
    double sunElevAngle     = 0.0;
    double sunAzimAngle     = 0.0;
    double lamda            = 1.5;
    int    iterNum          = 50;
    float  cellSize         = 1.0f;
    float  xllcorner        = 0.0f;
    float  yllcorner        = 0.0f;
    bool   haveSunElev      = false;
    bool   haveSunAzim      = false;
    bool   haveLamda        = false;
    bool   haveIterNum      = false;
    bool   initDemFlag      = false;

    if (argc == 1) usage();

    /* ---- parse arguments -------------------------------------- */
    int i = 1;
    while (i < argc)
    {
        if (argv[i][0] == '-')
        {
            if (strlen(argv[i]) < 2)
            {
                cerr << "Unknown flag '" << argv[i] << "'\n";
                usage();
            }
            char flag = argv[i][1];
            if (flag != 'h' && i + 1 >= argc)
            {
                cerr << "Flag -" << flag << " requires an argument.\n";
                usage();
            }
            switch (flag)
            {
            case 'e':
                sunElevAngle = atof(argv[++i]);
                haveSunElev  = true;
                break;
            case 'a':
                sunAzimAngle = atof(argv[++i]);
                haveSunAzim  = true;
                break;
            case 'l':
                lamda      = atof(argv[++i]);
                haveLamda  = true;
                break;
            case 'n':
                iterNum     = atoi(argv[++i]);
                haveIterNum = true;
                break;
            case 'd':
                strncpy(initDemName, argv[++i], sizeof(initDemName) - 1);
                initDemFlag = true;
                break;
            case 'c':
                cellSize = (float)atof(argv[++i]);
                break;
            case 'x':
                xllcorner = (float)atof(argv[++i]);
                break;
            case 'y':
                yllcorner = (float)atof(argv[++i]);
                break;
            case 'h':
                usage();
                break;
            default:
                cerr << "Unknown option -" << flag << "\n";
                usage();
            }
        }
        else
        {
            /* positional: first non-flag → input, second → output */
            if (inName[0] == '\0')
                strncpy(inName, argv[i], sizeof(inName) - 1);
            else if (outName[0] == '\0')
                strncpy(outName, argv[i], sizeof(outName) - 1);
            else
            {
                cerr << "Unexpected extra argument: " << argv[i] << "\n";
                usage();
            }
        }
        i++;
    }

    /* ---- validate required arguments -------------------------- */
    if (inName[0]  == '\0') { cerr << "Missing input TIFF filename.\n";  usage(); }
    if (outName[0] == '\0') { cerr << "Missing output TIFF filename.\n"; usage(); }
    if (!haveSunElev)  { cerr << "Missing -e (sun elevation angle).\n"; usage(); }
    if (!haveSunAzim)  { cerr << "Missing -a (sun azimuth  angle).\n";  usage(); }
    if (!haveLamda)    { cerr << "Missing -l (lambda).\n";              usage(); }
    if (!haveIterNum)  { cerr << "Missing -n (iteration count).\n";     usage(); }

    /* ---- echo parameters -------------------------------------- */
    cout << "************************************\n";
    cout << "Processing parameters:\n";
    cout << "  input          : " << inName        << "\n";
    cout << "  output         : " << outName       << "\n";
    cout << "  sunElevAngle   : " << sunElevAngle  << " deg\n";
    cout << "  sunAzimAngle   : " << sunAzimAngle  << " deg\n";
    cout << "  lambda         : " << lamda         << "\n";
    cout << "  iterNum        : " << iterNum       << "\n";
    cout << "  cellSize       : " << cellSize      << " m/px\n";
    cout << "  xllcorner      : " << xllcorner     << "\n";
    cout << "  yllcorner      : " << yllcorner     << "\n";
    if (initDemFlag)
        cout << "  initDem        : " << initDemName << "\n";
    else
        cout << "  initDem        : (flat zero surface)\n";
    cout << "************************************\n";

    /* ---- read input image ------------------------------------- */
    header inImgHead = {0};
    inImgHead.cellsize  = cellSize;
    inImgHead.xllcorner = xllcorner;
    inImgHead.yllcorner = yllcorner;

    float** inImg = readTIFF(inName, &inImgHead);
    int nrows = inImgHead.nrows;
    int ncols = inImgHead.ncols;

    cout << "Input image: " << ncols << " x " << nrows << " pixels\n";

    /* ---- read or initialise DEM ------------------------------- */
    float** initDem = fimage(nrows, ncols);
    if (!initDem) { cerr << "Cannot allocate initDem.\n"; exit(1); }

    if (initDemFlag)
    {
        header demHead = {0};
        float** rawDem = readTIFF(initDemName, &demHead);

        if (demHead.nrows != nrows || demHead.ncols != ncols)
        {
            cerr << "Initial DEM '" << initDemName
                 << "' size (" << demHead.ncols << "x" << demHead.nrows
                 << ") does not match image size ("
                 << ncols << "x" << nrows << ").\n";
            exit(1);
        }
        /* Copy luminance values (treated as normalised elevation) */
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                initDem[r][c] = rawDem[r][c];

        free(rawDem[0]);
        free(rawDem);
    }
    else
    {
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                initDem[r][c] = 0.0f;
    }

    /* ---- shape-from-shading ----------------------------------- */
    cout << "Running shape-from-shading (mini_sfs)...\n";
    mini_sfs(inImg, initDem,
             nrows, ncols,
             (double)cellSize,
             sunElevAngle, sunAzimAngle,
             lamda, iterNum);
    /* After mini_sfs returns, initDem holds the computed elevation. */

    /* ---- write TIFF output ------------------------------------ */
    header outHead = inImgHead;   /* inherit spatial metadata */
    writeTIFF(outName, initDem, &outHead);

    /* ---- clean up --------------------------------------------- */
    free(inImg[0]);
    free(inImg);
    free(initDem[0]);
    free(initDem);

    cout << "\nShape-from-shading complete.  Bye!\n";
    return 0;
}
