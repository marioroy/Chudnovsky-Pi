#line 2 "../src/util.h"
/* Input and output functions for GMP/MPIR supporting large data.

 * Copyright 2018 by Mario Roy (marioeroy at gmail dot com)

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

#ifndef UTIL_H
#define UTIL_H

#if defined(USE_GMP) || !defined(USE_MPIR)
 #include <gmp.h>
#else
 #include <mpir.h>
#endif

#define __ABS(x) __GMP_ABS(x)

/* Input, output for mpf_t missing in {gmp,mpir}.h.
 *
 * Return the number of bytes read or written. If an error occurs, or the
 * end-of-file is reached, the return value is a short count (or zero).
 */

size_t mpf_inp_raw (mpf_ptr x, FILE *f)
{
  mp_size_t nrecs, x_prec, x_size; mp_exp_t x_exp;
  size_t bytes = 0;

  nrecs  = fread(&x_prec, sizeof(mp_size_t), 1, f);
  nrecs += fread(&x_size, sizeof(mp_size_t), 1, f);
  bytes += sizeof(mp_size_t) * nrecs;

  nrecs  = fread(&x_exp, sizeof(mp_exp_t), 1, f);
  bytes += sizeof(mp_exp_t) * nrecs;

  mpz_t z; mpz_init(z);

  if (__ABS(x_size) > 0) {
    if (z->_mp_alloc != __ABS(x_size))
      _mpz_realloc(z, __ABS(x_size));

    nrecs  = fread(z->_mp_d, sizeof(mp_limb_t), __ABS(x_size), f);
    bytes += sizeof(mp_limb_t) * nrecs;

    z->_mp_size = x_size;
  }

  mpf_set_z(x, z);
  x->_mp_prec = x_prec;
  x->_mp_exp  = x_exp;
  mpz_clear(z);

  return bytes;
}

size_t mpf_out_raw (FILE *f, mpf_srcptr x)
{
  mp_size_t nrecs, x_prec, x_size; mp_exp_t x_exp;
  size_t bytes = 0;

  x_prec = x->_mp_prec;
  x_size = x->_mp_size;
  x_exp  = x->_mp_exp;

  nrecs  = fwrite(&x_prec, sizeof(mp_size_t), 1, f);
  nrecs += fwrite(&x_size, sizeof(mp_size_t), 1, f);
  bytes += sizeof(mp_size_t) * nrecs;

  nrecs  = fwrite(&x_exp, sizeof(mp_exp_t), 1, f);
  bytes += sizeof(mp_exp_t) * nrecs;

  if (__ABS(x_size) > 0) {
    nrecs  = fwrite(x->_mp_d, sizeof(mp_limb_t), __ABS(x_size), f);
    bytes += sizeof(mp_limb_t) * nrecs;
  }

  return bytes;
}

/* Input, output for mpz_t, 8 bytes versus 4 for size.
 *   See https://gmplib.org/manual/Raw-Output-Internals.html
 *
 * Return the number of bytes read or written. If an error occurs, or the
 * end-of-file is reached, the return value is a short count (or zero).
 */

size_t mpz_inp_raw (mpz_ptr x, FILE *f)
{
  mp_size_t nrecs, x_size, x_alloc;
  size_t bytes = 0;

  nrecs  = fread(&x_size,  sizeof(mp_size_t), 1, f);
  nrecs += fread(&x_alloc, sizeof(mp_size_t), 1, f);
  bytes += sizeof(mp_size_t) * nrecs;

  if (x->_mp_alloc != x_alloc)
    _mpz_realloc(x, x_alloc);

  if (x_alloc > 0) {
    nrecs  = fread(x->_mp_d, sizeof(mp_limb_t), x_alloc, f);
    bytes += sizeof(mp_limb_t) * nrecs;
  }

  x->_mp_size = x_size;

  return bytes;
}

size_t mpz_out_raw (FILE *f, mpz_srcptr x)
{
  mp_size_t nrecs, x_size, x_alloc;
  size_t bytes = 0;

  x_size  = x->_mp_size;
  x_alloc = x->_mp_alloc;

  nrecs  = fwrite(&x_size,  sizeof(mp_size_t), 1, f);
  nrecs += fwrite(&x_alloc, sizeof(mp_size_t), 1, f);
  bytes += sizeof(mp_size_t) * nrecs;

  if (x_alloc > 0) {
    nrecs  = fwrite(x->_mp_d, sizeof(mp_limb_t), x_alloc, f);
    bytes += sizeof(mp_limb_t) * nrecs;
  }

  return bytes;
}

#endif /* UTIL_H */

