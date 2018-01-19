#line 2 "../src/perl-chudnovsky.c"
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

#include "pgmp-chudnovsky.h"
#include "util.h"

#if defined(_WIN32) && (defined(__GNUC__) && __GNUC_VERSION__ < 40800)
  typedef unsigned long pthread_t;
#else
# include <pthread.h>
#endif

mpz_t p0, q0, g0;

typedef struct {
  double cputime;
  mpz_t *x1, *x2;
} thr_mul_t;

void *_mul (void *thr_arg)
{
  double t = wall_clock();
  thr_mul_t *thr_data = (thr_mul_t *) thr_arg;

  mpz_t *x1 = thr_data->x1;
  mpz_t *x2 = thr_data->x2;

  mpz_mul(*x1, *x1, *x2);

  thr_data->cputime = wall_clock()-t;

  return ((void *) 0);
}

uint_t chudnovsky_init (SV *digits_sv)
{
  uint64_t digits;

  #ifdef __LP64__
  digits = SvUV(digits_sv);
  #else
  digits = strtoull(SvPV_nolen(digits_sv), NULL, 10);
  #endif

  /* Turn off safe signal handling so that C code responds immediately
   * when pressing Ctrl-C or receiving a signal.
   */
  PL_signals |= PERL_SIGNALS_UNSAFE_FLAG;

  uint_t terms = digits / DIGITS_PER_ITER;
  mpf_set_default_prec((long)(digits * BITS_PER_DIGIT + 16));

  mpz_init(p0); mpz_init(q0); mpz_init(g0);

  return terms;
}

uint64_t chudnovsky_max_digits ()
{
  return MAX_DIGITS;
}

void chudnovsky_build_sieve (uint_t terms)
{
  char *args[] = { NULL };
  call_argv("sieve_begin_time", G_DISCARD|G_VOID, args);

  sieve_size = max(3*5*23*29+1, terms*6);
  sieve = (sieve_t *) malloc(sizeof(sieve_t)*sieve_size/2);
  build_sieve(sieve_size, sieve);

  call_argv("sieve_end_time", G_DISCARD|G_VOID, args);
}

void chudnovsky_free_sieve ()
{
  free(sieve);
}

void chudnovsky_bs (uint_t a, uint_t b, uint_t terms, uint_t level, int fd_i, uint_t depth)
{
  fac_t fp1, fg1, ftmp, fmul;
  mpz_t gcd, p1, q1, g1;
  uint_t j;

  fac_init(fp1), fac_init(ftmp), mpz_init(gcd);
  fac_init(fg1), fac_init(fmul);

  tmp_t *tmp = (tmp_t *) malloc(sizeof(tmp_t) * (depth - 1));

  for (j = 0; j < depth - 1; j++) {
    mpz_init(tmp[j].p),  mpz_init(tmp[j].q),  mpz_init(tmp[j].g);
    fac_init(tmp[j].fp), fac_init(tmp[j].fg), tmp[j].cleared = 0;
  }

  if (a == 0) {
    bs(p0, q0, g0, fp1, fg1, a, b, terms, level, gcd, ftmp, fmul, tmp, 0, 1);
  } else {
    mpz_init(p1), mpz_init(q1), mpz_init(g1);
    bs(p1, q1, g1, fp1, fg1, a, b, terms, level, gcd, ftmp, fmul, tmp, 0, 1);
  }

  for (j = 0; j < depth - 1; j++) {
    if (!tmp[j].cleared) {
      mpz_clear(tmp[j].p),  mpz_clear(tmp[j].q),  mpz_clear(tmp[j].g);
      fac_clear(tmp[j].fp), fac_clear(tmp[j].fg);
    }
  }

  free(tmp);

  fac_clear(fp1), fac_clear(ftmp), mpz_clear(gcd);
  fac_clear(fg1), fac_clear(fmul);

  if (a > 0) {
    FILE *file_i = fdopen(fd_i, "wb");

    mpz_out_raw(file_i, p1); mpz_clear(p1);
    mpz_out_raw(file_i, q1); mpz_clear(q1);

    if (b < terms)
      mpz_out_raw(file_i, g1);

    mpz_clear(g1);

    fflush(file_i);
    fclose(file_i);
  }
}

double chudnovsky_sum (uint_t i, uint_t k, int fd_i, int fd_k, char *path_k, int gflag)
{
  double join_begin, pthread_time = 0.0;
  mpz_t p1, q1, g1, p2, q2, g2;
  mpz_t *p, *q, *g;

  mpz_init(p1), mpz_init(q1), mpz_init(g1);
  mpz_init(p2), mpz_init(q2), mpz_init(g2);

  /* load k */

  FILE *file_k = fdopen(fd_k, "r+b");
  mpz_inp_raw(p2, file_k);
  mpz_inp_raw(q2, file_k);

  if (gflag)
    mpz_inp_raw(g2, file_k);

  fclose(file_k);
  unlink(path_k);

  /* load i */

  FILE *file_i = fdopen(fd_i, "r+b");

  if (i == 0) {
    /* inside data thr/proc where data resides */
    p = &p0, q = &q0, g = &g0;
  } else {
    /* other, read from disk, save to disk later */
    p = &p1, q = &q1, g = &g1;
    mpz_inp_raw(p1, file_i);
    mpz_inp_raw(q1, file_i);
    mpz_inp_raw(g1, file_i);
  }

  rewind(file_i); (void) ftruncate(fd_i, 0);
  fflush(file_i);

  /* sum */

  pthread_t thr1 = 0, thr2 = 0;
  thr_mul_t thr1_mul, thr2_mul;

  thr1_mul.x1 = p;    // mpz_mul(*p, *p, p2)
  thr1_mul.x2 = &p2;
  thr1_mul.cputime = 0.0;

  thr2_mul.x1 = q;    // mpz_mul(*q, *q, p2)
  thr2_mul.x2 = &p2;
  thr2_mul.cputime = 0.0;

#if defined(_WIN32) && (defined(__GNUC__) && __GNUC_VERSION__ < 40800)
  /* On the Windows platform, run serially if compiled using older GCC */
  _mul((void *) &thr1_mul);
  _mul((void *) &thr2_mul);

#else
  /* If reached ulimit -u threshold, run serially silently */
  if (pthread_create(&thr1, NULL, _mul, (void *) &thr1_mul))
    thr1 = 0, _mul((void *) &thr1_mul);
  if (pthread_create(&thr2, NULL, _mul, (void *) &thr2_mul))
    thr2 = 0, _mul((void *) &thr2_mul);

#endif

  mpz_mul(q2, q2, *g);

  join_begin = wall_clock();

  if (thr2 && !pthread_join(thr2, NULL))
    pthread_time += thr2_mul.cputime;
  if (thr1 && !pthread_join(thr1, NULL))
    pthread_time += thr1_mul.cputime;

  pthread_time -= wall_clock() - join_begin;

  mpz_clear(p2);
  mpz_add(*q, *q, q2);
  mpz_clear(q2);

  if (gflag)
    mpz_mul(*g, *g, g2);
  else
    mpz_clear(*g);

  mpz_clear(g2);

  /* save */

  if (i > 0) mpz_out_raw(file_i, p1);
  mpz_clear(p1);

  if (i > 0) mpz_out_raw(file_i, q1);
  mpz_clear(q1);

  if (gflag) {
    if (i > 0) mpz_out_raw(file_i, g1);
    mpz_clear(g1);
  }

  fflush(file_i);
  fclose(file_i);

  return pthread_time;
}

void chudnovsky_sqrt (int fd_c)
{
  mpf_t ci; mpf_init(ci);
  FILE *file_c = fdopen(fd_c, "wb");

  my_sqrt_ui(ci, C);
  mpf_out_raw(file_c, ci);
  mpf_clear(ci);

  fflush(file_c);
  fclose(file_c);
}

void chudnovsky_final (uint64_t digits, int out, uint_t terms, char *path_c)
{
  mpf_t pi, qi, ci;
  uint64_t psize, qsize;
  time_t now; struct tm *localtm;

  if (terms == 0) {
    mpz_set_ui(p0, 1);
    mpz_set_ui(q0, 0);
  }

  /*
     prepare to convert integers to floats

           p*(C/D)*sqrt(C)
     pi = -----------------
              (q+A*p)
  */
  psize = mpz_sizeinbase(p0, 10);
  qsize = mpz_sizeinbase(q0, 10);

  mpz_addmul_ui(q0, p0, A);
  mpz_mul_ui(p0, p0, C/D);

  mpf_init(pi), mpf_set_z(pi, p0);
  mpz_clear(p0);

  mpf_init(qi), mpf_set_z(qi, q0);
  mpz_clear(q0);

  /* final step */

  char *args[] = { NULL };
  call_argv("div_begin_time", G_DISCARD|G_VOID, args);

  my_div(qi, pi, qi);
  mpf_clear(pi);

  mpf_init(ci);
  call_argv("wait_sqrt",      G_DISCARD|G_VOID, args);
  FILE *file_c = fopen(path_c, "rb");
  mpf_inp_raw(ci, file_c);

  fclose(file_c);
  unlink(path_c);

  call_argv("div_end_time",   G_DISCARD|G_VOID, args);
  call_argv("mul_begin_time", G_DISCARD|G_VOID, args);

  mpf_mul(qi, qi, ci);
  mpf_clear(ci);

  call_argv("mul_end_time",   G_DISCARD|G_VOID, args);
  call_argv("total_time",     G_DISCARD|G_VOID, args);

  fprintf(stderr,
    "# P size = %llu digits (%f)\n"
    "# Q size = %llu digits (%f)\n",
    (unsigned long long) psize, (double) psize/digits,
    (unsigned long long) qsize, (double) qsize/digits);

  now = time(0), localtm = localtime(&now);

  fprintf(stderr,"#   end date = %s\n", asctime(localtm));
  fflush(stderr);

  /* output Pi */

  if (out == 1) {
    mp_exp_t exp = 0; int i;
    char *str = malloc(digits+16+2);

    mpf_get_str(&(str[1]), &exp, 10, digits+16, qi);
    for (i = 0; i < exp; i++) str[i] = str[i+1];

    str[exp] = '.', str[digits+2] = 0;
    fwrite(str, sizeof(char), digits+2, stdout), fflush(stdout);
    fprintf(stderr, "\n"), fflush(stderr);

    free((void *) str);
  }
  else if (out >= 2 && out <= 14) {
    output_digits(qi, digits, out);
  }

  mpf_clear(qi);
}

