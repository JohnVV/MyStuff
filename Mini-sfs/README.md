

- to build 

g++ -O2 -Ddebug=0 mini.cpp mini_sfs.cpp fft2d.cpp imageio.cpp -I. -lm -ltiff -lfftw3f_threads -lfftw3f -lpthread -o Mini_sfs


------------------------------------------------------

Usage:
  mini [options] <input.tif> <output.tif>

Required options:

  -e <deg>   Sun elevation angle (decimal degrees)

  -a <deg>   Sun azimuth  angle  (decimal degrees, N=0 CW)
  
  -l <val>   Lambda – regularisation parameter (0.4 – 50)
  
  -n <int>   Maximum iteration count (e.g. 50)

Optional:

  -d <file>  Initial DEM as TIFF (default: flat zero surface)
  
  -c <val>   Cell size in metres per pixel (default 1.0)
  
  -x <val>   Easting  of lower-left corner (default 0.0)
  
  -y <val>   Northing of lower-left corner (default 0.0)
  

-------------------------------------------------------
Example:

./Mini_sfs -e 35.0 -a 90.0 -l 1.5 -n 50 PIA21750.tiff PIA21750.DEM.tiff

--------------------------------------------------------
original code from this 2002/3 paper:

"Derivation of surface topography and terrain parameters from single satellite image using shape-from-shading technique"

https://www.sciencedirect.com/science/article/abs/pii/S0098300403001389#preview-section-introduction


updated to the current gcc

updated to use multi CPU cores  

updated the fft code to use the current fftw3_boost 
























