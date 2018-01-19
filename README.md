Chudnovsky-Pi
=============

A parallel demonstration for computing Pi digits in C and Perl using the fast [Chudnovsky](https://en.wikipedia.org/wiki/Chudnovsky_algorithm) algorithm. C code embeds [OpenMP](https://computing.llnl.gov/tutorials/openMP/) directives to consume many cores. Perl scripts `pi-hobo.pl` and `pi-thrs.pl` use [MCE::Hobo](https://metacpan.org/pod/MCE::Hobo) and [threads](https://metacpan.org/pod/threads) respectively, and [pthreads](https://computing.llnl.gov/tutorials/pthreads/) for 2nd-level workers.

# Content

```text
   bin/
     pi-gmp.exe          C executable using OpenMP and GMP
     pi-hobo.pl          Perl script using Inline::C and MCE::Hobo
     pi-mpir.exe         C executable using OpenMP and MPIR
     pi-thrs.pl          Perl script using Inline::C and threads

   lib/
     Perl modules        Inline, Inline::C, MCE, MCE::Shared, and
                         Parse::RecDescent
   src/
     Makefile            For building pi-gmp.exe and pi-mpir.exe
     extra/              Parallel recursion support for mpn_get_str
     perl-chudnovsky.c   Code used by Perl via Inline::C
     pgmp-chudnovsky.c   Code containing main and OpenMP directives
     pgmp-chudnovsky.h   Common code for perl/pgmp-chudnovsky.c
     typemap             Typemap configuration used by Inline::C
     util.h              Inp & out functions supporting large data
                           E.g. mpf/mpz_inp_raw, mpf/mpz_out_raw
```

# Dependencies

On Microsoft Windows, use [64-bit Cygwin](http://www.cygwin.com) to not have
limitations. Select `gcc-core`, `gcc-g++`, `libcrypt-devel`, `m4`, `make`,
`perl`, and `yasm` during installation.

* A C compiler, gmake or make, and m4. Yasm is required to build MPIR.
* [GMP 5.1.0a](https://gmplib.org) or later.
* [MPIR 2.6.0](http://mpir.org) or later.
* [YASM 1.2.0](http://yasm.tortall.net/Download.html) or later.

The following modules ship with Perl typically.

* [ExtUtils::MakeMaker](https://metacpan.org/pod/ExtUtils::MakeMaker)
* [Storable](https://metacpan.org/pod/Storable)
* [Time::HiRes](https://metacpan.org/pod/Time::HiRes)

# Building

* **GMP**

Strawberry Perl users may skip this step as GMP is included. By default,
`make install` will install files in `/usr/local`. You can specify an
installation prefix other than `/usr/local` using `--prefix`,
for instance `--prefix=$HOME`.

```text
   cd gmp-6.1.2
   ./configure --disable-static --enable-shared --enable-cxx
   make -j 4
   make check
   sudo make install
```

* **YASM**

The `yasm` package by the OS vendor is fine if 1.2.0 or later.

```text
   cd yasm-1.3.0
   ./configure
   make -j 4
   sudo make install
```

* **MPIR**

MPIR may provide better performance on Microsoft Windows.

```text
   cd mpir-3.0.0
   ./configure --disable-static --enable-shared --enable-cxx
   make -j 4
   make check
   sudo make install
```

* **CODE**

Library selection is possible via `CFLAGS`. The default builds against GMP.
Please use gmake when available.

```text
   CFLAGS="-O2 -DUSE_GMP"  (or)  CFLAGS="-O2 -DUSE_MPIR"
   Machines with 128 GiB: CFLAGS="-O2 -DUSE_GMP -DBIG_SIEVE=1"
     
   Strawberry Perl <  v5.26  dmake
   Strawberry Perl >= v5.26  gmake

   Cygwin    make  CFLAGS="-O2 -DUSE_MPIR"
   FreeBSD   gmake CC=clang CFLAGS="..."
   Linux     make
   Mac OS X  make
   Solaris   gmake note: src/extra supports GMP only on SPARC machines
   Other OS  make
```

Make without a target builds Inline C objects for Perl scripts residing in
`../bin/`. The C objects are stored in `../.Inline/`.

```text
   cd Chudnovsky-Pi-master/src

   make CFLAGS="-O2 -DBIG_SIEVE=0 -DUSE_GMP"    # default on all OS'es
   make CFLAGS="-O2 -DBIG_SIEVE=0 -DUSE_MPIR"   # preferred on Cygwin

   make pi-gmp    # builds the binary executable using GMP
   make pi-mpir   # builds the binary executable using MPIR
```

# Usage

The usage is similar for `pi-gmp.exe`, `pi-mpir.exe`, and `pi-thrs.pl`.

```text
   SYNOPSIS
       perl pi-hobo.pl <digits> [ <option> <threads> ]

       <digits>  digits of Pi to output

       <option>  0 - just run (default)
                 1 - output digits only
                 2 - output digits (2 columns)
                 3 - output digits (3 columns)
                 N - output digits (N columns, max 14)

       <threads> number of threads (default 1)
                 specify 'auto' to run on all cores

   EXAMPLES
       perl pi-hobo.pl 10000000 1 auto | md5sum
           bc3234ae2e3f6ec7737f037b375eabec  -

       perl pi-hobo.pl 100000000 1 auto | md5sum
           969bfe295b67da45b68086eb05a8b031  -

       perl pi-hobo.pl 100000000 5 auto > pi.txt
```

# Limitations

The following limitations apply to 32-bit OS'es and Strawberry Perl.

* Maximum 16 threads on Windows, excluding Cygwin
* 120 million digits for 32-bit binaries, all OS'es
* 640 million digits for 64-bit Strawberry Perl

To compute more than 640 million digits on Microsoft Windows, install
[64-bit Cygwin](http://www.cygwin.com). I tested 1 and 2 billion digits,
limited by available memory. Mutex locking using threads is slow in Cygwin
during mpz/mpf_init for temporary variables. Run `pi-hobo.pl` instead for
best performance.

```
   # on Cygwin
   perl pi-hobo.pl 1000000000 1 8 | md5sum
   3901670f41a84174103bd6a8f07651c0 *-

   perl pi-hobo.pl 2000000000 1 8 | md5sum
   dcf466792a8958becbb05b74b983d8b1 *-
```

# Acknowledgements

The code is derived from examples on the web.

* https://gmplib.org/list-archives/gmp-discuss/2008-November/003444.html
* https://gmplib.org/download/misc/gmp-chudnovsky.c

