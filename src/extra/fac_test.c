
/*
 * Output the factorial of n via mpz_out_str to standard output.
 *
 * Add -I/usr/local/include -L/usr/local/lib or other path if needed
 * macOS's native compiler lacks OpenMP support, use pthread instead
 *
 * Using OpenMP:
 *   gcc -DPARALLEL -DUSE_GMP  -O2 -fopenmp fac_test.c -lgmp  -lm -o fac_test
 *   gcc -DPARALLEL -DUSE_MPIR -O2 -fopenmp fac_test.c -lmpir -lm -o fac_test
 *
 * Using pthreads:
 *   gcc -DPARALLEL -DUSE_GMP  -O2 -pthread fac_test.c -lgmp  -lm -o fac_test
 *   gcc -DPARALLEL -DUSE_MPIR -O2 -pthread fac_test.c -lmpir -lm -o fac_test
 *
 * time ./fac_test 7200000 > out
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#if !defined(WIN32)
# include <sys/time.h>
#endif

#if defined(USE_GMP) || !defined(USE_MPIR)
# include <gmp.h>
# if defined(PARALLEL)
#  if defined(_OPENMP)
#   include "./gmp/mpn_get_str_omp.c"
#   include "./gmp/mpz_out_str.c"
#  else
#   include "./gmp/mpn_get_str_thr.c"
#   include "./gmp/mpz_out_str.c"
#  endif
# endif
#else
# include <mpir.h>
# if defined(PARALLEL)
#  if defined(_OPENMP)
#   include "./mpir/mpn_get_str_omp.c"
#   include "./mpir/mpz_out_str.c"
#  else
#   include "./mpir/mpn_get_str_thr.c"
#   include "./mpir/mpz_out_str.c"
#  endif
# endif
#endif


// clock_gettime isn't available on some platforms, e.g. Darwin
//
// https://blog.habets.se/2010/09/
//   gettimeofday-should-never-be-used-to-measure-time.html

#if defined(__GNUC__) && !defined(__GNUC_VERSION__)
# define __GNUC_VERSION__ (__GNUC__ * 10000 + __GNUC_MINOR__ * 100)
#endif

double wall_clock ()
{
#if !defined(CLOCK_MONOTONIC) || (defined(__GNUC__) && __GNUC_VERSION__ < 40800)
  struct timeval timeval;

  (void) gettimeofday (&timeval, (void *) 0);
  return (double) timeval.tv_sec +
         (double) timeval.tv_usec / 1000000.0;
#else
  struct timespec timeval;

  (void) clock_gettime (CLOCK_MONOTONIC, &timeval);
  return (double) timeval.tv_sec +
         (double) timeval.tv_nsec / 1000000000.0;
#endif
}


void fact (int n)
{
  double wbegin, wend;
  mpz_t  p;

  wbegin = wall_clock();
  mpz_fac_ui(p, n);
  wend = wall_clock();

  fprintf(stderr, "mpz_fac_ui  : %9.3f secs.\n", wend - wbegin);
  fflush(stderr);

  wbegin = wall_clock();
  mpz_out_str(stdout, 10, p);
  wend = wall_clock();

  fprintf(stderr, "mpz_out_str : %9.3f secs.\n", wend - wbegin);
  fflush(stderr);

  mpz_clear(p);
}


int main (int argc, char *argv[])
{
  int n;

  if (argc <= 1) {
    printf("Usage: %s <number>\n", argv[0]);
    return 1;
  }
  n = atoi(argv[1]);
  assert(n >= 0);
  fact(n);

  return 0;
}

