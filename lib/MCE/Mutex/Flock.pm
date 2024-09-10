###############################################################################
## ----------------------------------------------------------------------------
## MCE::Mutex::Flock - Mutex locking via Fcntl.
##
###############################################################################

package MCE::Mutex::Flock;

use strict;
use warnings;

no warnings qw( threads recursion uninitialized once );

our $VERSION = '1.900';

use base 'MCE::Mutex';
use Fcntl ':flock';
use Scalar::Util 'looks_like_number';
use Time::HiRes 'alarm';

my $tid = $INC{'threads.pm'} ? threads->tid() : 0;

sub CLONE {
    $tid = threads->tid() if $INC{'threads.pm'};
}

sub MCE::Mutex::Flock::_guard::DESTROY {
    my ($pid, $obj) = @{ $_[0] };
    CORE::flock ($obj->{_fh}, LOCK_UN), $obj->{ $pid } = 0 if $obj->{ $pid };

    return;
}

sub DESTROY {
    my ($pid, $obj) = ($tid ? $$ .'.'. $tid : $$, @_);
    $obj->unlock(), close(delete $obj->{_fh}) if $obj->{ $pid };
    unlink $obj->{path} if ($obj->{_init} && $obj->{_init} eq $pid);

    return;
}

sub _open {
    my ($pid, $obj) = ($tid ? $$ .'.'. $tid : $$, @_);
    return if exists $obj->{ $pid };

    open $obj->{_fh}, '+>>:raw:stdio', $obj->{path}
        or Carp::croak("Could not create temp file $obj->{path}: $!");

    return;
}

###############################################################################
## ----------------------------------------------------------------------------
## Public methods.
##
###############################################################################

my ($id, $prog_name) = (0);

$prog_name =  $0;
$prog_name =~ s{^.*[\\/]}{}g;
$prog_name =  'perl' if ($prog_name eq '-e' || $prog_name eq '-');

sub new {
    my ($class, %obj) = (@_, impl => 'Flock');

    if (! defined $obj{path}) {
        my ($pid, $tmp_dir, $tmp_file) = ( abs($$) );

        if ($ENV{TEMP} && -d $ENV{TEMP} && -w _) {
            if ($^O =~ /mswin|mingw|msys|cygwin/i) {
                $tmp_dir  = $ENV{TEMP};
                $tmp_dir .= ($^O eq 'MSWin32') ? "\\Perl-MCE" : "/Perl-MCE";
                mkdir $tmp_dir unless (-d $tmp_dir);
            }
            else {
                $tmp_dir = $ENV{TEMP};
            }
        }
        elsif ($ENV{TMPDIR} && -d $ENV{TMPDIR} && -w _) {
            $tmp_dir = $ENV{TMPDIR};
        }
        elsif (-d '/tmp' && -w _) {
            $tmp_dir = '/tmp';
        }
        else {
            Carp::croak("No writable dir found for a temp file");
        }

        $id++, $tmp_dir =~ s{[\\/]$}{};

        # remove tainted'ness from $tmp_dir
        if ($^O eq 'MSWin32') {
            ($tmp_file) = "$tmp_dir\\$prog_name.$pid.$tid.$id" =~ /(.*)/;
        } else {
            ($tmp_file) = "$tmp_dir/$prog_name.$pid.$tid.$id" =~ /(.*)/;
        }

        $obj{_init} = $tid ? $$ .'.'. $tid : $$;
        $obj{ path} = $tmp_file.'.lock';

        # test open
        open my $fh, '+>>:raw:stdio', $obj{path}
            or Carp::croak("Could not create temp file $obj{path}: $!");

        close $fh;

        # set permission
        chmod 0600, $obj{path};
    }
    else {
        # test open
        open my $fh, '+>>:raw:stdio', $obj{path}
            or Carp::croak("Could not obtain flock on file $obj{path}: $!");

        close $fh;
    }

    return bless(\%obj, $class);
}

sub lock {
    my ($pid, $obj) = ($tid ? $$ .'.'. $tid : $$, shift);
    $obj->_open() unless exists $obj->{ $pid };

    CORE::flock ($obj->{_fh}, LOCK_EX), $obj->{ $pid } = 1
        unless $obj->{ $pid };

    return;
}

sub guard_lock {
    &lock(@_);
    bless([ $tid ? $$ .'.'. $tid : $$, $_[0] ], MCE::Mutex::Flock::_guard::);
}

*lock_exclusive = \&lock;

sub lock_shared {
    my ($pid, $obj) = ($tid ? $$ .'.'. $tid : $$, shift);
    $obj->_open() unless exists $obj->{ $pid };

    CORE::flock ($obj->{_fh}, LOCK_SH), $obj->{ $pid } = 1
        unless $obj->{ $pid };

    return;
}

sub unlock {
    my ($pid, $obj) = ($tid ? $$ .'.'. $tid : $$, shift);

    CORE::flock ($obj->{_fh}, LOCK_UN), $obj->{ $pid } = 0
        if $obj->{ $pid };

    return;
}

sub synchronize {
    my ($pid, $obj, $code) = ($tid ? $$ .'.'. $tid : $$, shift, shift);
    my (@ret);

    return unless ref($code) eq 'CODE';

    $obj->_open() unless exists $obj->{ $pid };

    # lock, run, unlock - inlined for performance
    my $guard = bless([ $pid, $obj ], MCE::Mutex::Flock::_guard::);
    unless ($obj->{ $pid }) {
        CORE::flock ($obj->{_fh}, LOCK_EX), $obj->{ $pid } = 1;
    }
    (defined wantarray)
      ? @ret = wantarray ? $code->(@_) : scalar $code->(@_)
      : $code->(@_);

    return wantarray ? @ret : $ret[-1];
}

*enter = \&synchronize;

sub timedwait {
    my ($obj, $timeout) = @_;
    die 'MCE::Mutex::Flock::timedwait() unimplemented in this platform'
        if ($^O eq 'MSWin32');

    $timeout = 1 unless defined $timeout;
    Carp::croak('MCE::Mutex: timedwait (timeout) is not valid')
        if (!looks_like_number($timeout) || $timeout < 0);

    $timeout = 0.0003 if $timeout < 0.0003;

    local $@; local $SIG{ALRM} = sub { alarm 0; die "timed out\n" };
    eval { alarm $timeout; $obj->lock_exclusive };
    alarm 0;

    ( $@ && $@ eq "timed out\n" ) ? '' : 1;
}

1;

__END__

###############################################################################
## ----------------------------------------------------------------------------
## Module usage.
##
###############################################################################

=head1 NAME

MCE::Mutex::Flock - Mutex locking via Fcntl

=head1 VERSION

This document describes MCE::Mutex::Flock version 1.900

=head1 DESCRIPTION

A Fcntl implementation for C<MCE::Mutex>.

The API is described in L<MCE::Mutex>.

=over 3

=item new

=item lock

=item lock_exclusive

=item lock_shared

=item guard_lock

=item unlock

=item synchronize

=item enter

=item timedwait

=back

=head1 AUTHOR

Mario E. Roy, S<E<lt>marioeroy AT gmail DOT comE<gt>>

=cut

