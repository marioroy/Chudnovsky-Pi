#!/usr/bin/env perl
###############################################################################
# -----------------------------------------------------------------------------
# Pi computation using Chudnovsky's algortithm.
#
# AUTHOR
#     By Mario Roy using Inline::C, MCE::Shared, and MCE::Hobo.
#     Thank you, Dana Jacobsen (@danaj) for help with GMP.
#
# SEE ALSO
#     https://metacpan.org/pod/Math::Prime::Util::GMP
#     https://github.com/davidcarver/Parallel-Chudnovsky-PI
#     https://gmplib.org/download/misc/gmp-chudnovsky.c
#     https://gmplib.org/pi-with-gmp.html
#
###############################################################################

use strict;
use warnings;
use 5.010001;

use Cwd qw/ abs_path /;
use Time::HiRes qw/ time /;

my ( $prog_name, $prog_dir, $base_dir, $CFLAGS, $MPLIB );

BEGIN {
   $prog_name = $0;            $prog_name =~ s{^.*[\\/]}{}g;
   $prog_dir  = abs_path($0);  $prog_dir  =~ s{[\\/][^\\/]*$}{};
   $base_dir  = $prog_dir;     $base_dir  =~ s{[\\/][^\\/]*$}{};

   no warnings 'uninitialized';

   if ( @ARGV == 0 || $ARGV[1] eq 'auto' ) {
      printf {*STDERR} "\n".
         "SYNOPSIS\n".
         "    perl $prog_name <digits> [ <option> <threads> ]\n".
         "\n".
         "    <digits>  digits of Pi to output\n".
         "\n".
         "    <option>  0 - just run (default)\n".
         "              1 - output digits only\n".
         "              2 - output digits (2 columns)\n".
         "              3 - output digits (3 columns)\n".
         "              N - output digits (N columns, max 14)\n".
         "\n".
         "    <threads> number of threads (default 1)\n".
         "              specify 'auto' to run on all cores\n".
         "\n".
         "EXAMPLES\n".
         "    perl $prog_name 10000000 1 auto | md5sum\n".
         "        bc3234ae2e3f6ec7737f037b375eabec  -\n".
         "\n".
         "    perl $prog_name 100000000 1 auto | md5sum\n".
         "        969bfe295b67da45b68086eb05a8b031  -\n".
         "\n".
         "    perl $prog_name 100000000 5 auto > pi.txt\n".
         "\n";

      exit 1;
   }

   # The CFLAGS and MPLIB vars are relevant the first time only.
   # To reset, build again by passing argument, CFLAGS="...".

   if ( @ARGV && index($ARGV[0], 'CFLAGS') == 0 ) {
      $CFLAGS = substr($ARGV[0], 7), require File::Path;
      File::Path::rmtree("${base_dir}/.Inline");
   } else {
      $CFLAGS = '-O2 -DBIG_SIEVE=0 -DUSE_GMP';
   }

   $MPLIB = ( index($CFLAGS,'-DUSE_GMP' ) >= 0 ||
              index($CFLAGS,'-DUSE_MPIR') <  0 ) ? 'gmp' : 'mpir';

   mkdir "${base_dir}/.Inline" unless -d "${base_dir}/.Inline";
   $ENV{PERL_INLINE_DIRECTORY} = "${base_dir}/.Inline";

   unshift @INC, "$base_dir/lib";
}

###############################################################################
# -----------------------------------------------------------------------------
# Chudnovsky Pi.
#
###############################################################################

package c;

# C functions to be built under the c namespace.

# Inline::C is loaded here and inside a BEGIN block so to run
# before MCE::Signal's import subroutine. Any errors reported
# by the C compiler will cause the script to terminate, before
# making the tmp_dir.

BEGIN {
   use Inline "C" => Config =>
      TYPEMAPS => "${base_dir}/src/typemap",
      CCFLAGSEX => "${CFLAGS} -Wno-attributes -fno-stack-protector",
      inc => "-I${base_dir}/src -I/usr/local/include",
      libs => "-L/usr/local/lib -l${MPLIB} -lm",
      clean_after_build => 1;

   use Inline "C" => "${base_dir}/src/perl-chudnovsky.c";

   if ( @ARGV && index($ARGV[0],'CFLAGS') == 0 ) {
      exit 0;  # called from Makefile
   }
}

package main;

# Make a temporary dir under /dev/shm if writable, otherwise
# under /tmp. On the Windows platform, make the dir in the
# user's home direcory:
#    c:\Users\<userid>\AppData\Local\Temp\Perl-MCE\<temp-dir>

use MCE::Signal qw/ $tmp_dir -use_dev_shm /;

use MCE::Hobo 1.831;   # 1.831 minimally
use MCE::Mutex;
use MCE::Shared;

my $digits  = shift // 100;
my $output  = shift // 0;
my $threads = shift // 1;

my ( $cputime, $total_cputime, $total_wallclock );
my ( $ncpus, $max_digits, $mutex, $dataq );

$ncpus      = MCE::Util::get_ncpu();
$threads    = $ncpus if $threads eq 'auto';
$max_digits = c::chudnovsky_max_digits();

if ( $^O eq 'MSWin32' && $ncpus > 16 ) {
   $ncpus = 16; # limit for Windows only
}

if ( $digits > $max_digits ) {
   print {*STDERR} "Number of digits reset from $digits to $max_digits\n";
   $digits = $max_digits;
}

chudnovsky_pi( $digits, $output, $threads );

exit 0;

sub chudnovsky_pi {
   my ( $digits, $output, $threads ) = @_;

   my $terms = c::chudnovsky_init($digits);
   my $depth = 0;
   my $begin;

   if ( $threads < 1 || ( $terms <= 0 && $threads > 1 ) ) {
      print {*STDERR} "Number of threads reset from $threads to 1\n";
      $threads = 1;
   }
   elsif ( $terms > 0 && $terms < $threads && $threads <= $ncpus ) {
      print {*STDERR} "Number of threads reset from $threads to $terms\n";
      $threads = $terms;
   }
   elsif ( $threads > $ncpus ) {
      print {*STDERR} "Number of threads reset from $threads to $ncpus\n";
      $threads = $ncpus;
   }

   # In the event IO::FDPass is not available, construct a shared-queue first
   # before constructing other shared-objects and/or calling MCE::Hobo->init.

   $dataq = MCE::Shared->queue( await => 1, fast => 1 );
   $mutex = MCE::Mutex->new();

   $total_cputime   = MCE::Shared->scalar( 0 );
   $total_wallclock = MCE::Shared->scalar( 0 );
   $cputime         = MCE::Shared->scalar( 0 );

   $depth++ while ( (1 << $depth) < $terms );
   $depth++;

   printf {*STDERR} "# start date = %s\n", scalar(localtime());
   printf {*STDERR} "# terms = %lu, depth = %lu, threads = %d, logical cores = %d\n",
      $terms, $depth, $threads, $ncpus;

   MCE::Hobo->init(
      max_workers => ($threads < $ncpus) ? $threads : $ncpus,
      posix_exit  => 1
   );

   # defined task handler

   my $task_handler = sub {
      my ( $task, @args ) = @_;

      if ( $task eq 'bs' ) {
         my ( $i, $a, $b, $terms, $cores_depth ) = @args;
         my ( $path_i, $wbegin ) = ( "$tmp_dir/$i", time() );

         open my $fh_i, ">:unix:raw", $path_i or die "error '$path_i': $!";
         c::chudnovsky_bs($a, $b, $terms, $cores_depth, fileno($fh_i), $depth);
         close $fh_i;

         $cputime->incrby(time() - $wbegin);
      }
      elsif ( $task eq 'sum' ) {
         my ( $i, $k, $gflag ) = @args;
         my ( $path_i, $wbegin ) = ( "$tmp_dir/".($i), time() );
         my ( $path_k ) = ( "$tmp_dir/".($i+$k) );

         open my $fh_i, "+<:unix:raw", $path_i or die "error '$path_i': $!";
         open my $fh_k, "+<:unix:raw", $path_k or die "error '$path_k': $!";

         my $pthread_time = c::chudnovsky_sum(
            $i, $k, fileno($fh_i), fileno($fh_k), $path_k, $gflag );

         close $fh_i; close $fh_k;

         $cputime->incrby((time() - $wbegin) + $pthread_time);
      }
      elsif ( $task eq 'sqrt' ) {
         my ( $path_c, $wbegin ) = ( "$tmp_dir/c", time() );

         open my $fh_c, ">:unix:raw", $path_c or die "error '$path_c': $!";
         c::chudnovsky_sqrt(fileno($fh_c));
         close $fh_c;

         $cputime->incrby(time() - $wbegin);
      }
      elsif ( $task eq 'final' ) {
         my ( $digits, $output, $terms ) = @args;
         my $path_c = "$tmp_dir/c";

         c::chudnovsky_final($digits, $output, $terms, $path_c);
      }
      elsif ( $task eq 'free' ) {
         c::chudnovsky_free_sieve();
      }

      return;
   };

   # data resides under data_thr where final will run from

   c::chudnovsky_build_sieve($terms);

   my $data_thr = MCE::Hobo->create( sub {
      while ( defined ( my $task = $dataq->dequeue() ) ) {
         $task_handler->(@{ $task });
      }
   });

   if ( $terms == 0 ) {
      display_time('bs', 0.0, 0.0);
      display_time('sum', 0.0, 0.0);
   }
   else {
      my ( $cores_depth, $cores_size ) = ( 0 );
      my $mid = int($terms / $threads);
      my @thrs;

      $cores_depth++ while ((1 << $cores_depth) < $threads);
      $cores_size = 2 ** $cores_depth;

      # binary split

      $cputime->set(0), $begin = time();

      for my $i ( 0 .. $threads - 1 ) {
         my ( $a, $b ) = ( $i < $threads - 1 )
            ? ( $i * $mid, ($i + 1) * $mid )
            : ( $i * $mid, $terms );

         ( $i > 0 && $ncpus > 1 )
            ? push @thrs, MCE::Hobo->create(
                 $task_handler, 'bs', $i, $a, $b, $terms, $cores_depth
              )
            : $dataq->enqueue([ 'bs', $i, $a, $b, $terms, $cores_depth ]);
      }

      shift( @thrs )->join() while @thrs;

      # Note: To prevent the OS from performing a copy-on-write,
      # free the sieve resource first, inside the data-process.
      # This applies to non-threads only.

      $dataq->enqueue([ 'free' ]) if $^O ne 'MSWin32';
      $dataq->enqueue([ 'noop' ]);   # append non-task item
      $dataq->await(0);              # wait for thr_data worker

      # free sieve resource, main-process where sieve originated

      c::chudnovsky_free_sieve();

      display_time('bs', $cputime->get(), time() - $begin);

      # sum

      $cputime->set(0), $begin = time();

      for ( my $k = 1; $k < $cores_size; $k *= 2 ) {
         for ( my $i = 0; $i < $threads; $i = $i+2*$k ) {
            next if ( $i+$k >= $threads );
            my $gflag = ( $i+2*$k < $threads ) ? 1 : 0;

            ( $i > 0 && $ncpus > 1 )
               ? push @thrs, MCE::Hobo->create(
                    $task_handler, 'sum', $i, $k, $gflag
                 )
               : $dataq->enqueue([ 'sum', $i, $k, $gflag ]);
         }

         shift( @thrs )->join() while @thrs;
      }

      $dataq->enqueue([ 'noop' ]);   # append non-task item
      $dataq->await(0);              # wait for thr_data

      ( $threads > 1 )
         ? display_time('sum', $cputime->get(), time() - $begin)
         : display_time('sum', 0.0, 0.0);

      unlink "$tmp_dir/0";
   }

   # final step

   my $sqrt_thr;

   $cputime->set(0); $mutex->lock();

   ( $ncpus > 1 && $threads > 1 )
      ? $sqrt_thr = MCE::Hobo->create($task_handler, 'sqrt')
      : $dataq->enqueue([ 'sqrt' ]);

   $dataq->enqueue([ 'final', $digits, $output, $terms ]);
   $sqrt_thr->join() if defined $sqrt_thr;

   $mutex->unlock();    # release the lock, sqrt completed
   $dataq->end();       # terminate the queue
   $data_thr->join();   # reap the data thread
}

sub wait_sqrt {
   $mutex->lock();      # called by data_thr, wait for sqrt_thr
   $mutex->unlock();
}

###############################################################################
# -----------------------------------------------------------------------------
# Time routines called from C code.
#
###############################################################################

my $begin_time = 0;

sub display_time {
   my ( $desc, $cputime, $wallclock ) = @_;
   $cputime = $wallclock if ( $cputime < $wallclock );

   if ( $cputime > 0.0 && $wallclock > 0.0 ) {
      printf {*STDERR}
         "  %-8s  cputime = %9.2fs  wallclock = %8.2fs  factor = %5.1f\n",
         $desc, $cputime, $wallclock, $cputime / $wallclock;

   } else {
      $cputime = $wallclock = ($cputime > 0.0) ? $cputime : 0.0;

      printf {*STDERR}
         "  %-8s  cputime = %9.2fs  wallclock = %8.2fs  factor = %5.1f\n",
         $desc, $cputime, $wallclock, 1.0;
   }

   if ( $desc eq 'total' ) {
      printf {*STDERR} "%21s %9.2fm %12s %8.2fm\n",
         "", $cputime / 60.0, "", $wallclock / 60.0;
   } else {
      $total_cputime->incrby($cputime);
      $total_wallclock->incrby($wallclock);
   }
}

sub div_begin_time   { $begin_time = time(); }
sub mul_begin_time   { $begin_time = time(); }
sub sieve_begin_time { $begin_time = time(); }

sub div_end_time {
   my $wallclock = time() - $begin_time;
   my $sqrt_time = $cputime->get();
   my $cputime   = $wallclock + $sqrt_time;

   $wallclock += $sqrt_time if ( $ncpus == 1 || $threads == 1 );
   display_time('div/sqrt', $cputime, $wallclock);
}

sub mul_end_time {
   my $wallclock = time() - $begin_time;
   display_time('mul', $wallclock, $wallclock);
}

sub sieve_end_time {
   my $wallclock = time() - $begin_time;
   display_time('sieve', $wallclock, $wallclock);
}

sub total_time {
   display_time('total', $total_cputime->get(), $total_wallclock->get());
}

