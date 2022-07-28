/*  prime5.cpp - print number of decimal digits in biggest known prime number
    By Anthony C. Hay, 2013. See http://howtowriteaprogram.blogspot.com/

    This is free and unencumbered public domain software; see http://unlicense.org/
    I believe this code to be correct, but I may be wrong; use at your own risk.

    Compiled under Linux/macOS with g++ 4.2 or later
    macOS's native compiler lacks OpenMP support, use pthread instead

      Using OpenMP for parallel recursion in mpn_get_str:

        g++ -O3 -DPARALLEL -DTEST5 -fopenmp prime_test.cpp \
            -I/usr/local/include -L/usr/local/lib \
            -l mpir -o prime5_test -Wno-attributes

      Using pthreads for parallel recursion in mpn_get_str:

        g++ -O3 -DPARALLEL -DTEST5 -pthread prime_test.cpp \
            -I/usr/local/include -L/usr/local/lib \
            -l mpir -o prime5_test -Wno-attributes

    Compiled on Windows with VC++ 2010 using MPIR

      cl /nologo /EHs /O2 /DPARALLEL /DTEST5 \
          /I C:\mpir-2.6.0\lib\Win32\Release \
          prime_test.cpp /Feprime5_test.exe \
          /link C:\mpir-2.6.0\lib\Win32\Release\mpir.lib
*/

#include <iostream>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <stdint.h>
#include <mpir.h>

#if defined(PARALLEL)
# if defined(_OPENMP)
#  include "./mpir/mpn_get_str_omp.c"
#  include "./mpir/mpf_get_str.c"
#  include "./mpir/mpz_get_str.c"
# else
#  include "./mpir/mpn_get_str_thr.c"
#  include "./mpir/mpf_get_str.c"
#  include "./mpir/mpz_get_str.c"
# endif
#endif


// our arbitrary length unsigned number will be represented by a vector of
// "fragments", each fragment will be one of these
typedef mp_limb_t num_frag_t;
#define NUM_FRAG_T_SIZE GMP_LIMB_BITS
const int num_frag_t_size = NUM_FRAG_T_SIZE;

// arbitrary length unsigned number; least significant bits in lowest fragment
typedef std::vector<num_frag_t> num_vec_t;


// return given 'num' as a decimal string
std::string to_string(const num_vec_t & num)
{
    // use mpn_get_str() to do the binary to decimal conversion; this is
    // straight-forward because if we ensure our num_frag_t is the same
    // size as GMP's mp_limb_t then our binary representation will match
    // GMP's (what we called a fragment is what GMP calls a limb)

    num_vec_t n(num); // make a copy because mpn_get_str() destroys its argument
    n.push_back(0); // mpn_get_str() needs one extra high limb

    // allocate a buffer big enough to hold the whole decimal string
    // (log(2)/log(10) = 0.30102999566398119521373889472449)
    const size_t buf_len = static_cast<size_t>(num.size() * GMP_LIMB_BITS * 0.30103) + 3;
    std::vector<unsigned char> buf(buf_len);

    size_t len = mpn_get_str(&buf[0], 10, &n[0], num.size());

    // skip any leading zeros
    unsigned char * start = &buf[0];
    while (len && *start == 0) {
        --len;
        ++start;
    }
    // convert 0-9 -> '0'-'9'
    for (size_t i = 0; i < len; ++i)
        start[i] += '0'; // (C++ standard guarantees '0'..'9' are consecutive)

    return len ? std::string(start, start + len) : std::string("0");
}


// set given 'p' to the value of (2^n)-1 for given 'n'; requires n > 0
void make_prime(int n, num_vec_t & p)
{
    // a binary representation of (2^n)-1 is just an n-bit number with
    // every bit set; e.g. if n=4, (2^n)-1=15 decimal, or 1111 binary
    const int len = (n + num_frag_t_size - 1) / num_frag_t_size;
    const num_frag_t ones = ~static_cast<num_frag_t>(0);
    p = num_vec_t(len, ones);
    if (n % num_frag_t_size) {
        // the bits in p[len - 1] are already correct if n is a whole multiple
        // of num_frag_t_size (also, ones>>num_frag_t_size is undefined in C++)
        p[len - 1] = ones >> (num_frag_t_size - n % num_frag_t_size);
    }
}


// return a decimal string representation of (2^n)-1 for the given 'n'
std::string prime_str(int n)
{
    num_vec_t p;
    make_prime(n, p);
    return to_string(p);
}



#ifndef PRIME_UNDER_TEST

int main()
{
    // output the number of decimal digits in (2^57,885,161)-1, i.e. 17,425,170
    std::cout << prime_str(57885161).size() << '\n';
}

#endif
