/*------------------------------------------------------------------
  fft2d.cpp  –  2-D complex-to-complex FFT using FFTW 3.3
                Single-precision (fftwf_*) to match the float* buffers
                used throughout mini_sfs.cpp.
                Multi-core execution via libfftw3f_threads.

  Public interface (unchanged from the original):

    void fft2d(float* data, int nrows, int ncols, int isign);

      data   – row-major array of nrows*ncols interleaved (real,imag)
               pairs, i.e. length = nrows * ncols * 2.  The layout is
               identical to fftwf_complex* so no copy is needed.
      nrows  – number of rows    (any positive integer; FFTW3 is most
                                  efficient when the value is a product
                                  of small primes: 2, 3, 5, 7, 11, 13)
      ncols  – number of columns (same note as nrows)
      isign  – +1 : "forward"  in the original Numerical Recipes sense,
                    i.e. exp(+2πi·k/N)  → FFTW_BACKWARD
               -1 : "inverse"  in the original Numerical Recipes sense,
                    i.e. exp(−2πi·k/N)  → FFTW_FORWARD

  Normalisation: FFTW (like the original code) does NOT normalise.
  The caller (mini_sfs.cpp) divides by (nrows*ncols) explicitly after
  each inverse transform, exactly as before.

  Thread count (in decreasing priority):
    1. Environment variable MINI_FFT_THREADS  (e.g. export MINI_FFT_THREADS=4)
    2. std::thread::hardware_concurrency()    (all logical CPU cores)
    3. Fallback: 1 (single-threaded)

  Build:
    g++ ... -lfftw3f_threads -lfftw3f -lpthread -lm

  Install on Debian/Ubuntu:  sudo apt install libfftw3-dev
  Install on Fedora/RHEL  :  sudo dnf install fftw-devel
------------------------------------------------------------------*/

#include "mini.h"
#include <fftw3.h>
#include <thread>       /* std::thread::hardware_concurrency() */
#include <cstdlib>      /* getenv(), atoi()                    */

/* ----------------------------------------------------------------
   Thread initialisation
   ---------------------
   fftwf_init_threads() must be called exactly once before any other
   fftwf_* function.  We use a static local bool inside a helper so
   the call is guaranteed to happen first and only once regardless of
   the order of static initialisers across translation units.
---------------------------------------------------------------- */
namespace {

/* Returns the number of threads to use, resolving the three sources
   described in the header comment above. */
static int resolve_nthreads()
{
    /* 1. Environment variable override */
    const char* env = getenv("MINI_FFT_THREADS");
    if (env && atoi(env) > 0)
    {
        int n = atoi(env);
        cout << "fft2d: using " << n
             << " thread(s) (from MINI_FFT_THREADS)\n";
        return n;
    }

    /* 2. Hardware concurrency */
    unsigned int hw = std::thread::hardware_concurrency();
    int n = (hw > 0) ? static_cast<int>(hw) : 1;
    cout << "fft2d: using " << n << " thread(s) (hardware concurrency)\n";
    return n;
}

/* Called once at startup.  Returns the thread count actually set. */
static int init_fftw_threads()
{
    if (fftwf_init_threads() == 0)
    {
        fprintf(stderr, "fft2d: fftwf_init_threads() failed – "
                        "falling back to single-threaded execution.\n");
        return 1;
    }
    return resolve_nthreads();
}

/* Resolved once, at the time the first translation-unit static is
   initialised (i.e. before main() runs). */
static const int g_nthreads = init_fftw_threads();

/* ----------------------------------------------------------------
   Plan cache
   ----------
   fftwf_plan_with_nthreads(n) is a global setting that applies to
   all plans created after the call.  We call it once inside rebuild()
   immediately before fftwf_plan_dft_2d so the setting is always in
   effect when we need it.

   Plans are cached by dimension pair.  In this program the dimensions
   never change, so rebuild() is called exactly once.  If they ever do
   change, the old plans are destroyed and new ones created.
---------------------------------------------------------------- */
struct PlanCache
{
    int        nrows = 0;
    int        ncols = 0;
    fftwf_plan fwd   = nullptr;   /* FFTW_BACKWARD, isign = +1 */
    fftwf_plan inv   = nullptr;   /* FFTW_FORWARD,  isign = -1 */

    void rebuild(int nr, int nc)
    {
        destroy();

        /* Tell FFTW how many threads the next plan(s) should use.
           This must be called before every fftwf_plan_dft_2d call;
           it is a global, not per-plan, setting. */
        fftwf_plan_with_nthreads(g_nthreads);

        fftwf_complex* tmp =
            fftwf_alloc_complex(static_cast<size_t>(nr * nc));
        if (!tmp) {
            fprintf(stderr, "fft2d: fftwf_alloc_complex failed "
                            "(%d x %d).\n", nr, nc);
            exit(1);
        }

        fwd = fftwf_plan_dft_2d(nr, nc, tmp, tmp,
                                 FFTW_BACKWARD, FFTW_ESTIMATE);
        inv = fftwf_plan_dft_2d(nr, nc, tmp, tmp,
                                 FFTW_FORWARD,  FFTW_ESTIMATE);
        fftwf_free(tmp);

        if (!fwd || !inv) {
            fprintf(stderr, "fft2d: fftwf_plan_dft_2d failed "
                            "(%d x %d).\n", nr, nc);
            exit(1);
        }

        nrows = nr;
        ncols = nc;
    }

    void destroy()
    {
        if (fwd) { fftwf_destroy_plan(fwd); fwd = nullptr; }
        if (inv) { fftwf_destroy_plan(inv); inv = nullptr; }
        nrows = ncols = 0;
    }

    ~PlanCache()
    {
        destroy();
        /* Release all thread-local resources allocated by FFTW. */
        fftwf_cleanup_threads();
    }
};

static PlanCache g_cache;

} /* anonymous namespace */

/* ----------------------------------------------------------------- */
void fft2d(float* data, int nrows, int ncols, int isign)
{
    /* Rebuild plans only when dimensions change (normally once). */
    if (nrows != g_cache.nrows || ncols != g_cache.ncols)
        g_cache.rebuild(nrows, ncols);

    /*
      The float* buffer is laid out as consecutive (real, imag) pairs,
      which is exactly the memory layout of fftwf_complex (float[2]).
      A reinterpret_cast is therefore both correct and zero-copy.
    */
    fftwf_complex* cdata = reinterpret_cast<fftwf_complex*>(data);

    /*
      Sign convention mapping
      -----------------------
      Original  isign = +1  =>  exp(+2pi*i*k/N)  =>  FFTW_BACKWARD
      Original  isign = -1  =>  exp(-2pi*i*k/N)  =>  FFTW_FORWARD
    */
    fftwf_plan plan = (isign == 1) ? g_cache.fwd : g_cache.inv;

    /*
      fftwf_execute_dft executes the plan on a new data pointer.
      This is the recommended FFTW3 pattern for reusable plans
      (FFTW manual section 4.6: "New-array Execute Functions").
      When the plan was created with fftwf_plan_with_nthreads(n > 1),
      execute automatically spawns the worker threads internally;
      no change to the call site is needed.
    */
    fftwf_execute_dft(plan, cdata, cdata);
}
