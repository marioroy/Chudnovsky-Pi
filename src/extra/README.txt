//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Motivation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parallel recursion for mpn_get_str supporting OpenMP and pthreads.
By Mario Roy, 2018. See https://github.com/marioroy/binary-to-decimal.

Dear GMP/MPIR developers,

A common wish on the web is for mpn_get_str to run faster. Please, feel
free to disregard my humble attempt. For really "big" numbers, it still
takes a long time before reaching the initial divide-and-conquer inside
mpn_dc_get_str. At which point 2 threads run, then 4, and not to exceed
8 threads max.

Acknowledgement
  https://github.com/anthay/binary-to-decimal, by Anthony Hay

  prime_test.cpp is useful for validating changes in mpn/get_str.c
  added prime5.cpp (using MPIR) and prime6.cpp (using GMP)

Best,
  Mario Roy

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Enable parallel by including relevant header files inside C/C++ code.
// Likewise, it requires the appropiate CFLAGS option.
//
//   gcc/g++ -fopenmp  (or)  gcc/g++ -pthread
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// GMP: tested with gmp-5.1.0a minimally through gmp-6.1.2
//      also, benefits mpfr using gmp 5.1.0a or later

#include <gmp.h>

#if defined(_OPENMP)
# include "extra/gmp/mpn_get_str_omp.c"
#else
# include "extra/gmp/mpn_get_str_thr.c"
#endif

#include "extra/gmp/mpf_get_str.c"
#include "extra/gmp/mpz_get_str.c"
#include "extra/gmp/mpf_out_str.c"
#include "extra/gmp/mpz_out_str.c"

// MPIR: tested with mpir-2.6.0 minimally through mpir-3.0.0

#include <mpir.h>

#if defined(_OPENMP)
# include "extra/mpir/mpn_get_str_omp.c"
#else
# include "extra/mpir/mpn_get_str_thr.c"
#endif

#include "extra/mpir/mpf_get_str.c"
#include "extra/mpir/mpz_get_str.c"
#include "extra/mpir/mpf_out_str.c"
#include "extra/mpir/mpz_out_str.c"

// One may choose mpn_get_str_thr.c instead for OpenMP. That works too.

#if defined(_OPENMP)
# include "extra/{gmp,mpir}/mpn_get_str_thr.c"
#else
# ...
#endif

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ls -R extra/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

README.txt  gmp  mpir

extra/gmp:
COPYING         mpf_get_str.c      mpn_get_str_thr.c
gmp-impl.h      mpf_out_str.c      mpz_get_str.c
longlong.h      mpn_get_str_omp.c  mpz_out_str.c

extra/mpir:
COPYING         longlong.h         mpn_get_str_omp.c  mpz_out_str.c
arm             mpf_get_str.c      mpn_get_str_thr.c  x86
gmp-impl.h      mpf_out_str.c      mpz_get_str.c      x86_64

extra/mpir/arm:
longlong.h

extra/mpir/x86:
longlong.h

extra/mpir/x86_64:
longlong.h

The unchanged files longlong.h and mpf/mpz {get,out}_str.c are taken from
baseline. They are placed here so that building and runtime pick up mpf/mpz
{get,out}_str and subsequently mpn get_str from here versus lib{gmp,mpir}.

The gmp-impl.h file is greatly shrunked down. Only the relevant bits remain
for building successfully.

If desired, updating the baseline tree requires just one file change.
  cp extra/gmp/mpn_get_str_thr.c gmp-6.1.2/mpn/get_str.c
  cp extra/mpir/mpn_get_str_thr.c mpir-3.0.0/mpn/get_str.c

That requires adding -pthread to CFLAGS somewhere. The pthread solution
runs on platforms including languages lacking OpenMP support. E.g. Perl.

Thank you for GMP/MPIR.

// Cheers.

