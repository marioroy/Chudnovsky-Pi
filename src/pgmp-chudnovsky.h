#line 2 "../src/pgmp-chudnovsky.h"
/* Pi computation using Chudnovsky's algorithm.

 * Copyright 2002, 2005 Hanhong Xue (macroxue at yahoo dot com)

 * Slightly modified 2005 by Torbjorn Granlund (tege at swox dot com) to allow
   more than 2G digits to be computed.

 * Modified 2008 by David Carver (dcarver at tacc dot utexas dot edu) to enable
   multi-threading using the algorithm from "Computation of High-Precision
   Mathematical Constants in a Combined Cluster and Grid Environment" by
   Daisuke Takahashi, Mitsuhisa Sato, and Taisuke Boku.

 * Modified 2018 by Mario Roy (marioeroy at gmail dot com) to allow sieve to
   run with half memory consumption. Enabled terms and gflag optimizations
   for "bs" and "sum" respectively. Enabled 2nd-level parallelization in sum.
   Enabled parallel recursion in mpn_dc_get_str, see extra/README.txt.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
 * EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PGMP_CHUDNOVSKY_H
#define PGMP_CHUDNOVSKY_H

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <unistd.h>
#include <time.h>

#if !defined(_WIN32)
#include <sys/time.h>
#endif

#if defined(_WIN32)
#define strtoull _strtoui64
#endif

#if defined(USE_GMP) || !defined(USE_MPIR)
 #include <gmp.h>

 #if !defined(_WIN32)
  #if defined(_OPENMP)
   #include "extra/gmp/mpn_get_str_omp.c"
  #else
   #include "extra/gmp/mpn_get_str_thr.c"
  #endif
  #include "extra/gmp/mpf_get_str.c"
 #endif

#else
 #include <mpir.h>

 #if !defined(_WIN32)
  #if defined(_OPENMP)
   #include "extra/mpir/mpn_get_str_omp.c"
  #else
   #include "extra/mpir/mpn_get_str_thr.c"
  #endif
  #include "extra/mpir/mpf_get_str.c"
 #endif

#endif

#define BITS_PER_DIGIT   3.32192809488736234787  // log2(10)
#define DIGITS_PER_ITER  14.1816474627254776555  // log(53360^3)/log(10)
#define DOUBLE_PREC      53

#define A  13591409
#define B  545140134
#define C  640320
#define D  12

////////////////////////////////////////////////////////////////////////////

// clock_gettime isn't available on some platforms, e.g. Darwin
//
// https://blog.habets.se/2010/09/
//   gettimeofday-should-never-be-used-to-measure-time.html

#if defined(__GNUC__) && !defined(__GNUC_VERSION__)
#define __GNUC_VERSION__ (__GNUC__ * 10000 + __GNUC_MINOR__ * 100)
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

// For formatting comma-separated numbers, using the locale's thousands
// separator, if available. See http://c-faq.com/stdio/commaprint.html.

#define __MAXDIGITS (sizeof(int64_t) * 8 * sizeof(char) / 3) + 2

char *commify (uint64_t n)
{
  static int comma = '\0';
  static char retbuf[__MAXDIGITS];
  char *p = &retbuf[sizeof(retbuf)-1];
  int  i = 0;

  if (comma == '\0') {
    struct lconv *lcp = localeconv();
    if (lcp != NULL) {
      if (lcp->thousands_sep != NULL && *lcp->thousands_sep != '\0')
        comma = *lcp->thousands_sep;
      else
        comma = ',';
    }
  }

  *p = '\0';

  do {
    if (i%3 == 0 && i != 0) *--p = comma;
    *--p = '0' + n % 10;
    n /= 10, i++;
  } while (n != 0);

  return p;
}

// Display digits to standard output with spacing.

void output_digits (mpf_t pi, uint64_t digits, int columns)
{
  if (columns < 1) return;

  mp_exp_t exp = 0;
  uint64_t acc = 0;
  char *p, *str = malloc(digits+16+2);
  char *b, *buf = malloc(columns*11+1);
  int  acc_width, i, j, k, flag = 0, max = columns*10;

  snprintf(buf, __MAXDIGITS-1, "%s", commify(digits));
  acc_width = strlen(buf);

  mpf_get_str(&(str[1]), &exp, 10, digits+16, pi);
  for (i = 0; i < exp; i++) str[i] = str[i+1];

  str[exp] = '.', str[digits+2] = 0;
  b = buf, p = str, p += 2, i = j = 0;
  printf("3.");

  while (*p) {
    *b++ = *p++;

    if (++i % 10 == 0) {
      *b++ = ' ';

      if (i % max == 0) {
        *b = 0, acc += max;
        printf("%s :  %*s\n", buf, acc_width, commify(acc));
        if (++j % 10 == 0) { printf("\n"), j = 0; }
        if (*p) printf("  ");

        b = buf, flag = 1, i = 0;
      }
    }
  }

  if (i != 0 || digits == 0) {
    if (flag) {
      for (k = 10; k < max; k += 10) { if (i < k) *b++ = ' '; }
      *b = 0, acc += i;
      printf("%s %*s :  %*s\n", buf, max-i, "", acc_width, commify(acc));
    }
    else {
      if (i == 0 || i % 10 != 0) *b++ = ' ';
      *b = 0, acc += i;
      printf("%s :  %*s\n", buf, acc_width, commify(acc));
    }
  }

  fflush(stdout);

  if (!flag || j != 0)
    fprintf(stderr, "\n"), fflush(stderr);

  free((void *) buf);
  free((void *) str);
}

////////////////////////////////////////////////////////////////////////////

// sieve_t limit using uint64_t (-DBIG_SIEVE=1), default uint32_t
//
// GMP has a limit of 41 billion digits to not overflow mpz_t.
// The (uint32_t) sieve accomodates 10,151,618,680 digits max.
// On 32-bit HW, limit digits to consume less than 3.5 GiB.

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

#if __LP64__
 #if BIG_SIEVE
  static const uint64_t MAX_DIGITS = UINT64_C(40000000000);
  typedef uint64_t uint_t;
 #else
  static const uint64_t MAX_DIGITS = UINT64_C(10000000000);
  typedef uint32_t uint_t;
 #endif
#else
 #if defined(_WIN64)
  static const uint64_t MAX_DIGITS = UINT64_C(640000000);
  typedef uint32_t uint_t;
 #else
  static const uint64_t MAX_DIGITS = UINT64_C(120000000);
  typedef uint32_t uint_t;
 #endif
#endif

typedef struct {
  uint_t max_facs;
  uint_t num_facs;
  uint_t *fac;
  uint_t *pow;
} fac_t[1];

typedef struct {
  uint_t fac;
  uint_t pow;
  uint_t nxt;
} sieve_t;

sieve_t *sieve;
uint_t  sieve_size;

#define INIT_FACS 32

void fac_reset (fac_t f)
{
  f[0].num_facs = 0;
}

void fac_init_size (fac_t f, uint_t s)
{
  if (s < INIT_FACS)
    s = INIT_FACS;

  f[0].fac = malloc(s*sizeof(uint_t)*2);
  f[0].pow = f[0].fac + s;
  f[0].max_facs = s;

  fac_reset(f);
}

void fac_init (fac_t f)
{
  fac_init_size(f, INIT_FACS);
}

void fac_clear (fac_t f)
{
  free(f[0].fac);
}

void fac_resize (fac_t f, uint_t s)
{
  if (f[0].max_facs < s) {
    fac_clear(f);
    fac_init_size(f, s);
  }
}

/* f = base^pow */
void fac_set_bp (fac_t f, uint_t base, uint_t pow)
{
  uint_t i;

  for (i=0, base/=2; base>0; i++, base = sieve[base].nxt) {
    f[0].fac[i] = sieve[base].fac;
    f[0].pow[i] = sieve[base].pow*pow;
  }

  f[0].num_facs = i;
}

/* r = f*g */
void fac_mul2 (fac_t r, fac_t f, fac_t g)
{
  uint_t i, j, k;

  for (i=j=k=0; i<f[0].num_facs && j<g[0].num_facs; k++) {
    if (f[0].fac[i] == g[0].fac[j]) {
      r[0].fac[k] = f[0].fac[i];
      r[0].pow[k] = f[0].pow[i] + g[0].pow[j];
      i++; j++;
    } else if (f[0].fac[i] < g[0].fac[j]) {
      r[0].fac[k] = f[0].fac[i];
      r[0].pow[k] = f[0].pow[i];
      i++;
    } else {
      r[0].fac[k] = g[0].fac[j];
      r[0].pow[k] = g[0].pow[j];
      j++;
    }
  }
  for (; i<f[0].num_facs; i++, k++) {
    r[0].fac[k] = f[0].fac[i];
    r[0].pow[k] = f[0].pow[i];
  }
  for (; j<g[0].num_facs; j++, k++) {
    r[0].fac[k] = g[0].fac[j];
    r[0].pow[k] = g[0].pow[j];
  }

  r[0].num_facs = k;
}

/* f *= g */
void fac_mul (fac_t f, fac_t g, fac_t fmul)
{
  fac_t tmp;
  fac_resize(fmul, f[0].num_facs + g[0].num_facs);
  fac_mul2(fmul, f, g);
  tmp[0]  = f[0];
  f[0]    = fmul[0];
  fmul[0] = tmp[0];
}

/* f *= base^pow */
void fac_mul_bp (fac_t f, uint_t base, uint_t pow, fac_t ftmp, fac_t fmul)
{
  fac_set_bp(ftmp, base, pow);
  fac_mul(f, ftmp, fmul);
}

/* remove factors of power 0 */
void fac_compact (fac_t f)
{
  uint_t i, j;

  for (i=0, j=0; i<f[0].num_facs; i++) {
    if (f[0].pow[i]>0) {
      if (j<i) {
        f[0].fac[j] = f[0].fac[i];
        f[0].pow[j] = f[0].pow[i];
      }
      j++;
    }
  }

  f[0].num_facs = j;
}

/* convert factorized form to number */
void bs_mul (mpz_t r, uint_t a, uint_t b, fac_t fmul)
{
  uint_t i, j;

  if (b-a<=32) {
    mpz_set_ui(r, 1);
    for (i=a; i<b; i++)
      for (j=0; j<fmul[0].pow[i]; j++)
        mpz_mul_ui(r, r, fmul[0].fac[i]);
  } else {
    mpz_t r2;
    mpz_init(r2);
    bs_mul(r2, a, (a+b)/2, fmul);
    bs_mul(r, (a+b)/2, b, fmul);
    mpz_mul(r, r, r2);
    mpz_clear(r2);
  }
}

/* f /= gcd(f,g), g /= gcd(f,g) */
void fac_remove_gcd (mpz_t p, fac_t fp, mpz_t g, fac_t fg, mpz_t gcd, fac_t fmul)
{
  uint_t i, j, k, c;

  fac_resize(fmul, min(fp->num_facs, fg->num_facs));

  for (i=j=k=0; i<fp->num_facs && j<fg->num_facs; ) {
    if (fp->fac[i] == fg->fac[j]) {
      c = min(fp->pow[i], fg->pow[j]);
      fp->pow[i] -= c;
      fg->pow[j] -= c;
      fmul->fac[k] = fp->fac[i];
      fmul->pow[k] = c;
      i++; j++; k++;
    } else if (fp->fac[i] < fg->fac[j]) {
      i++;
    } else {
      j++;
    }
  }

  fmul->num_facs = k;

  if (fmul->num_facs) {
    bs_mul(gcd, 0, fmul->num_facs, fmul);
    mpz_divexact(p, p, gcd);
    mpz_divexact(g, g, gcd);
    fac_compact(fp);
    fac_compact(fg);
  }
}

void build_sieve(uint_t n, sieve_t *s)
{
  uint_t m, i, j, k;

  sieve_size = n;
  m = (uint_t) sqrt(n);
  memset(s, 0, sizeof(sieve_t)*n/2);

  s[1/2].fac = 1;
  s[1/2].pow = 1;

  for (i=3; i<=n; i+=2) {
    if (s[i/2].fac == 0) {
      s[i/2].fac = i;
      s[i/2].pow = 1;
      if (i<=m) {
        for (j=i*i, k=i/2; j<=n; j+=i+i, k++) {
          if (s[j/2].fac==0) {
            s[j/2].fac = i;
            if (s[k].fac == i) {
              s[j/2].pow = s[k].pow + 1;
              s[j/2].nxt = s[k].nxt;
            } else {
              s[j/2].pow = 1;
              s[j/2].nxt = k;
            }
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////

// r = sqrt(x)

void my_sqrt_ui (mpf_t r, uint64_t x)
{
  mpf_t t1, t2;
  uint64_t prec, bits, prec0;

  prec0 = mpf_get_prec(r);

  if (prec0<=DOUBLE_PREC) {
    mpf_set_d(r, sqrt(x));
    return;
  }

  bits = 0;
  for (prec=prec0; prec>DOUBLE_PREC;) {
    int bit = prec&1;
    prec = (prec+bit)/2;
    bits = bits*2+bit;
  }

  mpf_init(t1); mpf_init(t2);
  mpf_set_prec_raw(t1, DOUBLE_PREC);
  mpf_set_d(t1, 1/sqrt(x));

  while (prec<prec0) {
    prec *= 2;
    if (prec<prec0) {
      /* t1 = t1+t1*(1-x*t1*t1)/2; */
      mpf_set_prec_raw(t2, prec);
      mpf_mul(t2, t1, t1);         // half x half -> full
      mpf_mul_ui(t2, t2, x);
      mpf_ui_sub(t2, 1, t2);
      mpf_set_prec_raw(t2, prec/2);
      mpf_div_2exp(t2, t2, 1);
      mpf_mul(t2, t2, t1);         // half x half -> half
      mpf_set_prec_raw(t1, prec);
      mpf_add(t1, t1, t2);
    } else {
      break;
    }
    prec -= (bits&1);
    bits /= 2;
  }

  /* t2=x*t1, t1 = t2+t1*(x-t2*t2)/2; */
  mpf_set_prec_raw(t2, prec0/2);
  mpf_mul_ui(t2, t1, x);
  mpf_mul(r, t2, t2);              // half x half -> full
  mpf_ui_sub(r, x, r);
  mpf_mul(t1, t1, r);              // half x half -> half
  mpf_div_2exp(t1, t1, 1);
  mpf_add(r, t1, t2);

  mpf_clear(t1);
  mpf_clear(t2);
}

// r = y/x   WARNING: r cannot be the same as y.

#if __GNU_MP_RELEASE >= 50001 && !defined(__MPIR_RELEASE)
#define my_div mpf_div
#else
void my_div (mpf_t r, mpf_t y, mpf_t x)
{
  mpf_t t1, t2;
  uint64_t prec, bits, prec0;

  prec0 = mpf_get_prec(r);

  if (prec0<=DOUBLE_PREC) {
    mpf_set_d(r, mpf_get_d(y)/mpf_get_d(x));
    return;
  }

  bits = 0;
  for (prec=prec0; prec>DOUBLE_PREC;) {
    int bit = prec&1;
    prec = (prec+bit)/2;
    bits = bits*2+bit;
  }

  mpf_init(t1); mpf_init(t2);
  mpf_set_prec_raw(t1, DOUBLE_PREC);
  mpf_ui_div(t1, 1, x);

  while (prec<prec0) {
    prec *= 2;
    if (prec<prec0) {
      /* t1 = t1+t1*(1-x*t1); */
      mpf_set_prec_raw(t2, prec);
      mpf_mul(t2, x, t1);          // full x half -> full
      mpf_ui_sub(t2, 1, t2);
      mpf_set_prec_raw(t2, prec/2);
      mpf_mul(t2, t2, t1);         // half x half -> half
      mpf_set_prec_raw(t1, prec);
      mpf_add(t1, t1, t2);
    } else {
      prec = prec0;
      /* t2=y*t1, t1 = t2+t1*(y-x*t2); */
      mpf_set_prec_raw(t2, prec/2);
      mpf_mul(t2, t1, y);          // half x half -> half
      mpf_mul(r, x, t2);           // full x half -> full
      mpf_sub(r, y, r);
      mpf_mul(t1, t1, r);          // half x half -> half
      mpf_add(r, t1, t2);
      break;
    }
    prec -= (bits&1);
    bits /= 2;
  }

  mpf_clear(t1);
  mpf_clear(t2);
}
#endif

////////////////////////////////////////////////////////////////////////////

// binary splitting

typedef struct {
  mpz_t p, q, g;
  fac_t fp, fg;
  int cleared;
} tmp_t;

#define  pj tmp[j].p
#define  qj tmp[j].q
#define  gj tmp[j].g
#define fpj tmp[j].fp
#define fgj tmp[j].fg

void bs (mpz_t p1, mpz_t q1, mpz_t g1, fac_t fp1, fac_t fg1, uint_t a, uint_t b, uint_t terms, uint_t level, mpz_t gcd, fac_t ftmp, fac_t fmul, tmp_t *tmp, uint_t j, int clear_flag)
{
  uint_t i, mid;

  if (b-a == 1) {
    /*
      g(b-1,b) = (6b-5)(2b-1)(6b-1)
      p(b-1,b) = b^3 * C^3 / 24
      q(b-1,b) = (-1)^b*g(b-1,b)*(A+Bb).
    */
    mpz_set_ui(p1, b);
    mpz_mul_ui(p1, p1, b);
    mpz_mul_ui(p1, p1, b);
    mpz_mul_ui(p1, p1, (C/24) * (C/24));
    mpz_mul_ui(p1, p1, (C*24));

    mpz_set_ui(g1, 2*b-1);
    mpz_mul_ui(g1, g1, 6*b-1);
    mpz_mul_ui(g1, g1, 6*b-5);

    mpz_set_ui(q1, b);
    mpz_mul_ui(q1, q1, B);
    mpz_add_ui(q1, q1, A);
    mpz_mul   (q1, q1, g1);

    if (b % 2)
      mpz_neg(q1, q1);

    i = b; while ((i & 1) == 0) i >>= 1;

    fac_set_bp(fp1, i, 3);                  /* b^3 */
    fac_mul_bp(fp1, 3*5*23*29, 3, ftmp, fmul);
    fp1[0].pow[0]--;

    fac_set_bp(fg1, 2*b-1, 1);              /* 2b-1 */
    fac_mul_bp(fg1, 6*b-1, 1, ftmp, fmul);  /* 6b-1 */
    fac_mul_bp(fg1, 6*b-5, 1, ftmp, fmul);  /* 6b-5 */

  } else {
    /*
      p(a,b) = p(a,m) * p(m,b)
      g(a,b) = g(a,m) * g(m,b)
      q(a,b) = q(a,m) * p(m,b) + q(m,b) * g(a,m)
    */
    mid = a + (b-a) * 0.54;  /* tuning parameter */

    bs(p1, q1, g1, fp1, fg1, a, mid, terms, level+1,
      gcd, ftmp, fmul, tmp, j, 0);

    bs(pj, qj, gj, fpj, fgj, mid, b, terms, level+1,
      gcd, ftmp, fmul, tmp, j+1, 0);

    if (level >= 4)          /* tuning parameter */
      fac_remove_gcd(pj, fpj, g1, fg1, gcd, fmul);

    mpz_mul(qj, qj, g1);
    fac_mul(fp1, fpj, fmul);

    if (b < terms) {
      mpz_mul(g1, g1, gj);
      fac_mul(fg1, fgj, fmul);
    }
    if (clear_flag) {
      fac_clear(fpj), fac_clear(fgj), mpz_clear(gj);
      tmp[j].cleared = 1;
    }

    mpz_mul(q1, q1, pj);
    mpz_add(q1, q1, qj);

    if (clear_flag) mpz_clear(qj);

    mpz_mul(p1, p1, pj);

    if (clear_flag) mpz_clear(pj);
  }
}

#endif /* PGMP_CHUDNOVSKY_H */

