###############################################################################
## ----------------------------------------------------------------------------
## A hybrid LRU-plain cache helper class.
##
## An optimized, pure-Perl LRU implementation with extra performance when
## fetching items from the upper-section of the cache.
##
###############################################################################

package MCE::Shared::Cache;

use strict;
use warnings;

use 5.010001;

no warnings qw( threads recursion uninitialized numeric );

our $VERSION = '1.836';

## no critic (Subroutines::ProhibitExplicitReturnUndef)
## no critic (TestingAndDebugging::ProhibitNoStrict)

use Scalar::Util qw( dualvar );
use Time::HiRes qw( time );

use MCE::Shared::Base ();
use base 'MCE::Shared::Base::Common';
use bytes;

use constant {
   _DATA => 0,  # unordered data
   _KEYS => 1,  # LRU queue
   _INDX => 2,  # index into _KEYS
   _BEGI => 3,  # begin offset value
   _GCNT => 4,  # garbage count
   _EXPI => 5,  # max age, default disabled
   _SIZE => 6,  # max keys, default disabled
   _HREF => 7,  # for hash-like dereferencing
   _ITER => 8,  # for tied hash support
};

use overload (
   q("")    => \&MCE::Shared::Base::_stringify,
   q(0+)    => \&MCE::Shared::Base::_numify,
   q(%{})   => sub {
      no overloading;
      $_[0]->[_HREF] || do {
         # no circular reference to original, therefore no memory leaks
         tie my %h, __PACKAGE__.'::_href', bless([ @{ $_[0] } ], __PACKAGE__);
         $_[0]->[_HREF] = \%h;
      };
   },
   fallback => 1
);

no overloading;

###############################################################################
## ----------------------------------------------------------------------------
## TIEHASH, STORE, FETCH, DELETE, FIRSTKEY, NEXTKEY, EXISTS, CLEAR, SCALAR
##
###############################################################################

# TIEHASH ( max_keys => undef, max_age => undef );  # default
# TIEHASH ( { options }, @pairs );
# TIEHASH ( )

sub TIEHASH {
   my $class = shift;
   my $opts  = ( ref $_[0] eq 'HASH' ) ? shift : undef;

   if ( !defined $opts ) {
      $opts = {};
      for my $cnt ( 1 .. 3 ) {
         if ( @_ && $_[0] =~ /^(max_keys|max_age|_shared)$/ ) {
            $opts->{ $1 } = $_[1];
            splice @_, 0, 2;
         }
      }
   }

   my $size = MCE::Shared::Cache::_size( $opts->{'max_keys'} // undef );
   my $expi = MCE::Shared::Cache::_secs( $opts->{'max_age' } // undef );
   my ( $begi, $gcnt ) = ( 0, 0 );

   my $obj = bless [ {}, [], {}, \$begi, \$gcnt, \$expi, \$size ], $class;

   # The shared-server preallocates the two hashes automatically, where the
   # data resides, when constructed with MCE::Shared->cache( ... ).

   $obj->_prealloc() unless $opts->{_shared};
   $obj->mset(@_) if @_;

   $obj;
}

# STORE ( key, value )

sub STORE {
   my ( $data, $keys, $indx, $begi, $gcnt, $expi, $size ) = @{ $_[0] };

   # update existing key
   if ( defined ( my $off = $indx->{ $_[1] } ) ) {
      $off -= ${ $begi };

      # update expiration
      $keys->[ $off ] = dualvar( time + ${ $expi }, $_[1] )
         if ( defined ${ $expi } );

      # promote key if not last, inlined for performance
      if ( ! $off ) {
         return $data->{ $_[1] } = $_[2] if @{ $keys } == 1;

         push @{ $keys }, shift @{ $keys };
         $indx->{ $_[1] } = ++${ $begi } + @{ $keys } - 1;

         MCE::Shared::Cache::_gckeys_head( $keys, $begi, $gcnt )
            if ( ${ $gcnt } && !defined $keys->[0] );

         # safety to not overrun
         $_[0]->purge() if ( ${ $begi } > 2.147e9 );
      }
      elsif ( $off != @{ $keys } - 1 ) {
         push @{ $keys }, delete $keys->[ $off ];
         $indx->{ $_[1] } = ${ $begi } + @{ $keys } - 1;

         # GC keys if gcnt:size ratio is greater than 2:3
         $_[0]->purge() if ( ++${ $gcnt } > @{ $keys } * 0.667 );
      }

      return $data->{ $_[1] } = $_[2];
   }

   # insert key-value pair
   $data->{ $_[1] } = $_[2];
   $indx->{ $_[1] } = ${ $begi } + @{ $keys };

   push @{ $keys }, defined ${ $expi }
      ? dualvar( time + ${ $expi }, $_[1] )
      : "$_[1]";

   # evict the least used key, inlined for performance
   if ( defined ${ $size } && @{ $keys } - ${ $gcnt } > ${ $size } ) {
      my $key = shift @{ $keys };
      ${ $begi }++; delete $data->{ $key }; delete $indx->{ $key };

      MCE::Shared::Cache::_gckeys_head( $keys, $begi, $gcnt )
         if ( ${ $gcnt } && !defined $keys->[0] );

      # safety to not overrun
      $_[0]->purge() if ( ${ $begi } > 2.147e9 );
   }

   $_[2];
}

# FETCH ( key )

sub FETCH {

   # cache miss
   return undef if !defined ( my $off = $_[0]->[_INDX]{ $_[1] } );

   # cache hit
   my ( $data, $keys, $indx, $begi, $gcnt, $expi ) = @{ $_[0] };

   $off -= ${ $begi };

   # key expired
   $_[0]->del( $_[1] ), return undef if (
      defined ${ $expi } && $keys->[ $off ] < time
   );

   # promote key if not upper half, inlined for performance
   if ( ! $off ) {
      return $data->{ $_[1] } if @{ $keys } == 1;

      push @{ $keys }, shift @{ $keys };
      $indx->{ $_[1] } = ++${ $begi } + @{ $keys } - 1;

      MCE::Shared::Cache::_gckeys_head( $keys, $begi, $gcnt )
         if ( ${ $gcnt } && !defined $keys->[0] );

      # safety to not overrun
      $_[0]->purge() if ( ${ $begi } > 2.147e9 );
   }
   elsif ( $off - ${ $gcnt } < ((@{ $keys } - ${ $gcnt }) >> 1) ) {
      push @{ $keys }, delete $keys->[ $off ];
      $indx->{ $_[1] } = ${ $begi } + @{ $keys } - 1;

      # GC keys if gcnt:size ratio is greater than 2:3
      $_[0]->purge() if ( ++${ $gcnt } > @{ $keys } * 0.667 );
   }

   $data->{ $_[1] };
}

# DELETE ( key )

sub DELETE {
   my ( $data, $keys, $indx, $begi, $gcnt ) = @{ $_[0] };

   return undef if !defined ( my $off = delete $indx->{ $_[1] } );

   $off -= ${ $begi };

   # check the first key
   if ( ! $off ) {
      ${ $begi }++; shift @{ $keys };

      if ( ${ $gcnt } && !defined $keys->[0] ) {
         MCE::Shared::Cache::_gckeys_head( $keys, $begi, $gcnt );
      } elsif ( ! @{ $keys } ) {
         ${ $begi } = 0;
      }

      return delete $data->{ $_[1] };
   }

   # check the last key
   elsif ( $off == @{ $keys } - 1 ) {
      pop @{ $keys };

      if ( ${ $gcnt } && !defined $keys->[-1] ) {
         MCE::Shared::Cache::_gckeys_tail( $keys, $gcnt );
      } elsif ( ! @{ $keys } ) {
         ${ $begi } = 0;
      }

      return delete $data->{ $_[1] };
   }

   # must be a key somewhere in-between
   $keys->[ $off ] = undef;   # tombstone

   # GC keys if gcnt:size ratio is greater than 2:3
   $_[0]->purge() if ( ++${ $gcnt } > @{ $keys } * 0.667 );

   delete $data->{ $_[1] };
}

# FIRSTKEY ( )

sub FIRSTKEY {
   my $self = shift;
   $self->[_ITER] = [ $self->keys ];

   $self->NEXTKEY;
}

# NEXTKEY ( )

sub NEXTKEY {
   shift @{ $_[0]->[_ITER] };
}

# EXISTS ( key )

sub EXISTS {
   my ( $self, $key ) = @_;
   return '' if !defined ( my $off = $self->[_INDX]{ $key } );

   $self->del( $key ), return '' if (
      defined ${ $self->[_EXPI] } &&
      $self->[_KEYS][ $off -= ${ $self->[_BEGI] } ] < time
   );

   1;
}

# CLEAR ( )

sub CLEAR {
   my ( $data, $keys, $indx, $begi, $gcnt ) = @{ $_[0] };

   %{ $data } = @{ $keys } = %{ $indx } = ();
   ${ $begi } = ${ $gcnt } = 0;

   delete $_[0]->[_ITER];

   return;
}

# SCALAR ( )

sub SCALAR {
   $_[0]->_prune_head() if defined ${ $_[0]->[_EXPI] };

   scalar keys %{ $_[0]->[_DATA] };
}

###############################################################################
## ----------------------------------------------------------------------------
## Internal routines for preserving dualvar KEYS data during freeze-thaw ops.
##
###############################################################################

## Storable freeze-thaw

sub STORABLE_freeze {
   my ( $self, $cloning ) = @_;
   return if $cloning;

   my @TIME;

   if ( defined ${ $self->[_EXPI] } ) {
      for my $key ( @{ $self->[_KEYS] } ) {
         push @TIME, defined $key ? 0 + $key : $key;
      }
   }

   return '', [ @{ $self }, \@TIME ];
}

sub STORABLE_thaw {
   my ( $self, $cloning, $serialized, $ret ) = @_;
   return if $cloning;

   my $TIME = pop @{ $ret };
   @{ $self } = @{ $ret };

   if ( defined ${ $self->[_EXPI] } ) {
      my ( $i, $keys ) = ( 0, $self->[_KEYS] );

      for my $time ( @{ $TIME } ) {
         $keys->[$i] = dualvar( $time, $keys->[$i] ) if defined ( $time );
         $i++;
      }
   }

   return;
}

## Sereal freeze-thaw

sub FREEZE {
   my ( $self ) = @_;
   my @TIME;

   if ( defined ${ $self->[_EXPI] } ) {
      for my $key ( @{ $self->[_KEYS] } ) {
         push @TIME, defined $key ? 0 + $key : $key;
      }
   }

   return [ @{ $self }, \@TIME ];
}

sub THAW {
   my ( $class, $serializer, $data ) = @_;
   my $TIME = pop @{ $data };
   my $self = $class->new();

   @{ $self } = @{ $data };

   if ( defined ${ $self->[_EXPI] } ) {
      my ( $i, $keys ) = ( 0, $self->[_KEYS] );

      for my $time ( @{ $TIME } ) {
         $keys->[$i] = dualvar( $time, $keys->[$i] ) if defined ( $time );
         $i++;
      }
   }

   return $self;
}

###############################################################################
## ----------------------------------------------------------------------------
## _gckeys_head, _gckeys_tail, _inskey, _prealloc, _prune_head, _secs, _size
##
###############################################################################

# Called by MCE::Shared::Server::_share

sub _shared_init {
   $_[0]->_prealloc();
}

# GC start of list

sub _gckeys_head {
   my ( $keys, $begi, $gcnt ) = @_;
   my $i = 1;

   $i++ until ( defined $keys->[$i] );
   ${ $begi } += $i, ${ $gcnt } -= $i;
   splice @{ $keys }, 0, $i;

   return;
}

# GC end of list

sub _gckeys_tail {
   my ( $keys, $gcnt ) = @_;
   my $i = $#{ $keys } - 1;

   $i-- until ( defined $keys->[$i] );
   ${ $gcnt } -= $#{ $keys } - $i;
   splice @{ $keys }, $i + 1;

   return;
}

# insert or promote key

sub _inskey {
   my ( $data, $keys, $indx, $begi, $gcnt, $expi, $size ) = @{ $_[0] };

   # update existing key
   if ( defined ( my $off = $indx->{ $_[1] } ) ) {
      $off -= ${ $begi };

      # unset value if expired
      $data->{ $_[1] } = undef
         if ( defined ${ $expi } && $keys->[ $off ] < time );

      # update expiration
      $keys->[ $off ] = dualvar( time + ${ $expi }, $_[1] )
         if ( defined ${ $expi } );

      # promote key if not last, inlined for performance
      if ( ! $off ) {
         return if @{ $keys } == 1;

         push @{ $keys }, shift @{ $keys };
         $indx->{ $_[1] } = ++${ $begi } + @{ $keys } - 1;

         MCE::Shared::Cache::_gckeys_head( $keys, $begi, $gcnt )
            if ( ${ $gcnt } && !defined $keys->[0] );

         # safety to not overrun
         $_[0]->purge() if ( ${ $begi } > 2.147e9 );
      }
      elsif ( $off != @{ $keys } - 1 ) {
         push @{ $keys }, delete $keys->[ $off ];
         $indx->{ $_[1] } = ${ $begi } + @{ $keys } - 1;

         # GC keys if gcnt:size ratio is greater than 2:3
         $_[0]->purge() if ( ++${ $gcnt } > @{ $keys } * 0.667 );
      }

      return;
   }

   # insert key
   $indx->{ $_[1] } = ${ $begi } + @{ $keys };

   push @{ $keys }, defined ${ $expi }
      ? dualvar( time + ${ $expi }, $_[1] )
      : "$_[1]";

   # evict the least used key, inlined for performance
   if ( defined ${ $size } && @{ $keys } - ${ $gcnt } > ${ $size } ) {
      my $key = shift @{ $keys };
      ${ $begi }++; delete $data->{ $key }; delete $indx->{ $key };

      MCE::Shared::Cache::_gckeys_head( $keys, $begi, $gcnt )
         if ( ${ $gcnt } && !defined $keys->[0] );

      # safety to not overrun
      $_[0]->purge() if ( ${ $begi } > 2.147e9 );
   }

   return;
}

# preallocate two hashes internally

sub _prealloc {
   my ( $self ) = @_;
   my $size = ${ $self->[_SIZE] } // 500;

   keys %{ $self->[_DATA] } = $size + 1;
   keys %{ $self->[_INDX] } = $size + 1;

   return;
}

# prune start of list

sub _prune_head {
   return unless defined ${ $_[0]->[_EXPI] };

   my ( $data, $keys, $indx, $begi, $gcnt ) = @{ $_[0] };
   my ( $i, $time ) = ( 0, time );

   for my $k ( @{ $keys } ) {
      $i++, ${ $gcnt }--, next unless ( defined $k );
      last if ( $keys->[$i] > $time );

      delete $data->{ $k };
      delete $indx->{ $k };
      $i++;
   }

   ${ $begi } += $i, splice @{ $keys }, 0, $i if $i;

   return;
}

# compute seconds

sub _secs {
   my ( $secs ) = @_;

   # seconds, minutes, hours, days, weeks
   my %secs = ( '' => 1, s => 1, m => 60, h => 3600, d => 86400, w => 604800 );

   $secs = undef if ( defined $secs && $secs eq 'never' );

   if ( defined $secs ) {
      $secs = 0 if ( $secs eq 'now' );
      $secs = 0.0001 if ( $secs ne '0' && $secs < 0.0001 );
      $secs = $1 * $secs{ lc($2) } if $secs =~ /^(\d*\.?\d*)\s*([smhdw]?)/i;
   }

   $secs;
}

# compute size

sub _size {
   my ( $size ) = @_;

   # Digital Information Sizes Calculator
   # http://dr-lex.be/info-stuff/bytecalc.html

   # kibiBytes (KiB), mebiBytes (MiB)
   my %size = ( '' => 1, k => 1024, m => 1048576 );

   $size = undef if ( defined $size && $size eq 'unlimited' );

   if ( defined $size ) {
      $size = 0 if $size < 0;
      $size = $1 * $size{ lc($2) } if $size =~ /^(\d*\.?\d*)\s*([km]?)/i;
      $size = int( $size + 0.5 );
   }

   $size;
}

###############################################################################
## ----------------------------------------------------------------------------
## _find, iterator, keys, pairs, values
##
###############################################################################

# _find ( { getkeys => 1 }, "query string" )
# _find ( { getvals => 1 }, "query string" )
# _find ( "query string" ) # pairs

sub _find {
   my $self   = shift;
   my $params = ref($_[0]) eq 'HASH' ? shift : {};
   my $query  = shift;

   MCE::Shared::Base::_find_hash( $self->[_DATA], $params, $query, $self );
}

# iterator ( key [, key, ... ] )
# iterator ( "query string" )
# iterator ( )

sub iterator {
   my ( $self, @keys ) = @_;
   my $data = $self->[_DATA];

   if ( ! @keys ) {
      @keys = $self->keys;
   }
   elsif ( @keys == 1 && $keys[0] =~ /^(?:key|val)[ ]+\S\S?[ ]+\S/ ) {
      @keys = $self->keys($keys[0]);
   }
   elsif ( defined ${ $self->[_EXPI] } ) {
      $self->_prune_head();
   }

   return sub {
      return unless @keys;
      my $key = shift @keys;
      return ( $key => $data->{ $key } );
   };
}

# keys ( key [, key, ... ] )
# keys ( "query string" )
# keys ( )

sub keys {
   my $self = shift;

   $self->_prune_head() if defined ${ $self->[_EXPI] };

   if ( @_ == 1 && $_[0] =~ /^(?:key|val)[ ]+\S\S?[ ]+\S/ ) {
      $self->_find( { getkeys => 1 }, @_ );
   }
   else {
      if ( wantarray ) {
         my $data = $self->[_DATA];
         @_ ? map { exists $data->{ $_ } ? $_ : undef } @_
            : $self->_keys;
      }
      else {
         scalar CORE::keys %{ $self->[_DATA] };
      }
   }
}

# _keys ( )

sub _keys {
   my $self = shift;

   if ( ${ $self->[_EXPI] } ) {
      map { ''. $_ } ${ $self->[_GCNT] }
         ? grep defined($_), reverse @{ $self->[_KEYS] }
         : reverse @{ $self->[_KEYS] };
   }
   else {
      ${ $self->[_GCNT] }
         ? grep defined($_), reverse @{ $self->[_KEYS] }
         : reverse @{ $self->[_KEYS] };
   }
}

# pairs ( key [, key, ... ] )
# pairs ( "query string" )
# pairs ( )

sub pairs {
   my $self = shift;

   $self->_prune_head() if defined ${ $self->[_EXPI] };

   if ( @_ == 1 && $_[0] =~ /^(?:key|val)[ ]+\S\S?[ ]+\S/ ) {
      $self->_find( @_ );
   }
   else {
      if ( wantarray ) {
         my $data = $self->[_DATA];
         @_ ? map { $_ => $data->{ $_ } } @_
            : map { $_ => $data->{ $_ } } $self->_keys;
      }
      else {
         scalar CORE::keys %{ $self->[_DATA] };
      }
   }
}

# values ( key [, key, ... ] )
# values ( "query string" )
# values ( )

sub values {
   my $self = shift;

   $self->_prune_head() if defined ${ $self->[_EXPI] };

   if ( @_ == 1 && $_[0] =~ /^(?:key|val)[ ]+\S\S?[ ]+\S/ ) {
      $self->_find( { getvals => 1 }, @_ );
   }
   else {
      if ( wantarray ) {
         @_ ? @{ $self->[_DATA] }{ @_ }
            : @{ $self->[_DATA] }{ $self->_keys };
      }
      else {
         scalar CORE::keys %{ $self->[_DATA] };
      }
   }
}

###############################################################################
## ----------------------------------------------------------------------------
## assign, max_age, max_keys, mdel, mexists, mget, mset, peek, purge
##
###############################################################################

# assign ( key, value [, key, value, ... ] )

sub assign {
   $_[0]->clear; shift()->mset(@_);
}

# max_age ( [ secs ] )

sub max_age {
   my ( $self, $secs ) = @_;
   my $keys = $self->[_KEYS];
   my $expi = $self->[_EXPI];

   $secs = MCE::Shared::Cache::_secs( $secs );

   if ( defined $secs ) {
      $self->purge();
      if ( defined ${ $expi } ) {
         my $off    = $secs - ${ $expi };
         @{ $keys } = map { dualvar( $_ + $off, $_.'' ) } @{ $keys };
         ${ $expi } = $secs, $self->_prune_head();
      }
      else {
         my $time   = time;
         @{ $keys } = map { dualvar( $time + $secs, $_ ) } @{ $keys };
         ${ $expi } = $secs, $self->_prune_head();
      }
   }
   elsif ( @_ == 2 ) {
      $self->purge();
      if ( defined ${ $expi } ) {
         @{ $keys } = map { $_.'' } @{ $keys };
         ${ $expi } = undef;
      }
   }

   if ( defined wantarray ) {
      defined ${ $expi } ? ${ $expi } > 0 ? ${ $expi } : 'now' : 'never';
   }
}

# max_keys ( [ size ] )

sub max_keys {
   my ( $self, $size ) = @_;

   $size = MCE::Shared::Cache::_size( $size );

   if ( defined $size ) {
      my ( $data, $keys, $indx, $begi, $gcnt ) = @{ $self };
      my $count = CORE::keys( %{ $data } ) - $size;

      # evict the least used keys
      while ( $count-- > 0 ) {
         my $key = shift @{ $keys };
         ${ $begi }++; delete $data->{ $key }; delete $indx->{ $key };

         MCE::Shared::Cache::_gckeys_head( $keys, $begi, $gcnt )
            if ( ${ $gcnt } && !defined $keys->[0] );

         # safety to not overrun
         $self->purge() if ( ${ $begi } > 2.147e9 );
      }

      ${ $self->[_SIZE] } = $size;
   }
   elsif ( @_ == 2 ) {
      ${ $self->[_SIZE] } = undef;
   }

   if ( defined wantarray ) {
      defined ${ $self->[_SIZE] } ? ${ $self->[_SIZE] } : 'unlimited';
   }
}

# mdel ( key [, key, ... ] )

sub mdel {
   my $self = shift;
   my $cnt  = 0;

   while ( @_ ) {
      my $key = shift;
      $cnt++, $self->del( $key ) if $self->exists( $key );
   }

   $cnt;
}

# mexists ( key [, key, ... ] )

sub mexists {
   my $self = shift;

   while ( @_ ) {
      return '' unless $self->exists( shift );
   }

   1;
}

# mget ( key [, key, ... ] )

sub mget {
   my $self = shift;

   @_ ? map { $self->get( $_ ) } @_ : ();
}

# mset ( key, value [, key, value, ... ] )

sub mset {
   my $self = shift;

   while ( @_ ) {
      $self->set( splice( @_, 0, 2 ) );
   }

   defined wantarray ? $self->SCALAR : ();
}

# peek ( key )

sub peek {
   $_[0]->[_DATA]{ $_[1] };
}

# purge ( )

sub purge {
   my ( $data, $keys, $indx, $begi, $gcnt, $expi ) = @{ $_[0] };
   my $i; $i = ${ $begi } = ${ $gcnt } = 0;

   # purge in-place for minimum memory consumption

   if ( defined ${ $expi } ) {
      my $time = time;
      for my $k ( @{ $keys } ) {
         if ( defined $k && $k < $time ) {
            delete $data->{ $k };
            delete $indx->{ $k };
            next;
         }
         $keys->[ $i ] = $k, $indx->{ $k } = $i++ if ( defined $k );
      }
   }
   else {
      for my $k ( @{ $keys } ) {
         $keys->[ $i ] = $k, $indx->{ $k } = $i++ if ( defined $k );
      }
   }

   splice @{ $keys }, $i;

   return;
}

###############################################################################
## ----------------------------------------------------------------------------
## Sugar API, mostly resembles http://redis.io/commands#string primitives.
##
###############################################################################

# append ( key, string )

sub append {
   $_[0]->_inskey( $_[1] );
   length( $_[0]->[_DATA]{ $_[1] } .= $_[2] // '' );
}

# decr ( key )

sub decr {
   $_[0]->_inskey( $_[1] );
   --$_[0]->[_DATA]{ $_[1] };
}

# decrby ( key, number )

sub decrby {
   $_[0]->_inskey( $_[1] );
   $_[0]->[_DATA]{ $_[1] } -= $_[2] || 0;
}

# incr ( key )

sub incr {
   $_[0]->_inskey( $_[1] );
   ++$_[0]->[_DATA]{ $_[1] };
}

# incrby ( key, number )

sub incrby {
   $_[0]->_inskey( $_[1] );
   $_[0]->[_DATA]{ $_[1] } += $_[2] || 0;
}

# getdecr ( key )

sub getdecr {
   $_[0]->_inskey( $_[1] );
   $_[0]->[_DATA]{ $_[1] }-- // 0;
}

# getincr ( key )

sub getincr {
   $_[0]->_inskey( $_[1] );
   $_[0]->[_DATA]{ $_[1] }++ // 0;
}

# getset ( key, value )

sub getset {
   $_[0]->_inskey( $_[1] );

   my $old = $_[0]->[_DATA]{ $_[1] };
   $_[0]->[_DATA]{ $_[1] } = $_[2];

   $old;
}

# len ( key )
# len ( )

sub len {
   $_[0]->_prune_head() if defined ${ $_[0]->[_EXPI] };

   ( defined $_[1] )
      ? length $_[0]->get( $_[1] )
      : scalar CORE::keys %{ $_[0]->[_DATA] };
}

{
   no strict 'refs';

   *{ __PACKAGE__.'::new'    } = \&TIEHASH;
   *{ __PACKAGE__.'::set'    } = \&STORE;
   *{ __PACKAGE__.'::get'    } = \&FETCH;
   *{ __PACKAGE__.'::delete' } = \&DELETE;
   *{ __PACKAGE__.'::exists' } = \&EXISTS;
   *{ __PACKAGE__.'::clear'  } = \&CLEAR;
   *{ __PACKAGE__.'::del'    } = \&delete;
   *{ __PACKAGE__.'::merge'  } = \&mset;
   *{ __PACKAGE__.'::vals'   } = \&values;
}

# For on-demand hash-like dereferencing.

package # hide from rpm
   MCE::Shared::Cache::_href;

sub TIEHASH { $_[1] }

1;

__END__

###############################################################################
## ----------------------------------------------------------------------------
## Module usage.
##
###############################################################################

=head1 NAME

MCE::Shared::Cache - A hybrid LRU-plain cache helper class

=head1 VERSION

This document describes MCE::Shared::Cache version 1.836

=head1 DESCRIPTION

A cache helper class for use as a standalone or managed by L<MCE::Shared>.

This module implements a least-recently used (LRU) cache with its origin based
on L<MCE::Shared::Ordhash>, for its performance and low-memory characteristics.
It is both a LRU and plain implementation. LRU logic is applied to new items
and subsequent updates. A fetch however, involves LRU reorder only if the item
is found in the lower section of the cache. This equates to extra performance
for the upper section as fetch behaves similarly to accessing a plain cache.
Upon reaching its size restriction, it prunes items from the bottom of the
cache.

The 50% LRU-mode (bottom section), 50% plain-mode (upper-section) applies to
fetches only.

=head1 SYNOPSIS

 # non-shared or local construction for use by a single process

 use MCE::Shared::Cache;

 my $ca;

 $ca = MCE::Shared::Cache->new(); # max_keys => undef, max_age => undef
 $ca = MCE::Shared::Cache->new( { max_keys => 500 }, @pairs );

 $ca = MCE::Shared::Cache->new( max_keys => "unlimited", max_age => "never" );
 $ca = MCE::Shared::Cache->new( max_keys => undef, max_age => undef ); # ditto
 $ca = MCE::Shared::Cache->new( max_keys => 500, max_age => "1 hour" );
 $ca = MCE::Shared::Cache->new( max_keys => "4 KiB" ); # 4*1024
 $ca = MCE::Shared::Cache->new( max_keys => "1 MiB" ); # 1*1024*1024

 $ca = MCE::Shared::Cache->new( max_age  => "43200 seconds" );
 $ca = MCE::Shared::Cache->new( max_age  => 43200 );   # ditto
 $ca = MCE::Shared::Cache->new( max_age  => "720 minutes" );
 $ca = MCE::Shared::Cache->new( max_age  => "12 hours" );
 $ca = MCE::Shared::Cache->new( max_age  => "0.5 days" );
 $ca = MCE::Shared::Cache->new( max_age  => "1 week" );
 $ca = MCE::Shared::Cache->new( max_age  => "never" ); # or undef
 $ca = MCE::Shared::Cache->new( max_age  => "now" );   # or 0

 # construction for sharing with other threads and processes

 use MCE::Shared;

 my $ca;

 $ca = MCE::Shared->cache(); # max_keys => undef, max_age => undef
 $ca = MCE::Shared->cache( { max_keys => 500 }, @pairs );

 $ca = MCE::Shared->cache( max_keys => "unlimited", max_age => "never" );
 $ca = MCE::Shared->cache( max_keys => undef, max_age => undef ); # ditto
 $ca = MCE::Shared->cache( max_keys => 500, max_age => "1 hour" );
 $ca = MCE::Shared->cache( max_keys => "4 KiB" ); # 4*1024
 $ca = MCE::Shared->cache( max_keys => "1 MiB" ); # 1*1024*1024

 $ca = MCE::Shared->cache( max_age  => "43200 seconds" );
 $ca = MCE::Shared->cache( max_age  => 43200 );   # ditto
 $ca = MCE::Shared->cache( max_age  => "720 minutes" );
 $ca = MCE::Shared->cache( max_age  => "12 hours" );
 $ca = MCE::Shared->cache( max_age  => "0.5 days" );
 $ca = MCE::Shared->cache( max_age  => "1 week" );
 $ca = MCE::Shared->cache( max_age  => "never" ); # or undef
 $ca = MCE::Shared->cache( max_age  => "now" );   # or 0

 # hash-like dereferencing

 my $val = $ca->{$key};
 $ca->{$key} = $val;

 %{$ca} = ();

 # OO interface

 if ( !defined ( $val = $ca->get("some_key") ) ) {
    $val = $ca->set( some_key => "some_value" );
 }

 $val   = $ca->set( $key, $val );
 $val   = $ca->get( $key );
 $val   = $ca->delete( $key );              # del is an alias for delete
 $bool  = $ca->exists( $key );
 void   = $ca->clear();
 $len   = $ca->len();                       # scalar keys %{ $ca }
 $len   = $ca->len( $key );                 # length $ca->{ $key }

 $iter  = $ca->iterator( @keys );           # ($key, $val) = $iter->()
 @keys  = $ca->keys( @keys );               # @keys is optional
 %pairs = $ca->pairs( @keys );
 @vals  = $ca->values( @keys );             # vals is an alias for values

 $len   = $ca->assign( $key/$val pairs );   # equivalent to ->clear, ->mset
 $cnt   = $ca->mdel( @keys );
 @vals  = $ca->mget( @keys );
 $bool  = $ca->mexists( @keys );            # true if all keys exists
 $len   = $ca->mset( $key/$val pairs );     # merge is an alias for mset

 # included, sugar methods without having to call set/get explicitly

 $len   = $ca->append( $key, $string );     #   $val .= $string
 $val   = $ca->decr( $key );                # --$val
 $val   = $ca->decrby( $key, $number );     #   $val -= $number
 $val   = $ca->getdecr( $key );             #   $val--
 $val   = $ca->getincr( $key );             #   $val++
 $val   = $ca->incr( $key );                # ++$val
 $val   = $ca->incrby( $key, $number );     #   $val += $number
 $old   = $ca->getset( $key, $new );        #   $o = $v, $v = $n, $o

 # pipeline, provides atomicity for shared objects, MCE::Shared v1.09+

 @vals  = $ca->pipeline(                    # ( "a_a", "b_b", "c_c" )
    [ "set", foo => "a_a" ],
    [ "set", bar => "b_b" ],
    [ "set", baz => "c_c" ],
    [ "mget", qw/ foo bar baz / ]
 );

For normal hash behavior, the TIE interface is supported.

 # non-shared or local construction for use by a single process

 use MCE::Shared::Cache;

 tie my %ca, "MCE::Shared::Cache", max_keys => undef, max_age => undef;
 tie my %ca, "MCE::Shared::Cache", max_keys => 500, max_age => "1 hour";
 tie my %ca, "MCE::Shared::Cache", { max_keys => 500 }, @pairs;

 # construction for sharing with other threads and processes
 # one option is needed minimally to know to use MCE::Shared::Cache

 use MCE::Shared;

 tie my %ca, "MCE::Shared", max_keys => undef, max_age => undef;
 tie my %ca, "MCE::Shared", max_keys => 500, max_age => "1 hour";
 tie my %ca, "MCE::Shared", { max_keys => 500 }, @pairs;

 # usage

 my $val;

 if ( !defined ( $val = $ca{some_key} ) ) {
    $val = $ca{some_key} = "some_value";
 }

 $ca{some_key} = 0;

 tied(%ca)->incrby("some_key", 20);
 tied(%ca)->incrby(some_key => 20);

=head1 SYNTAX for QUERY STRING

Several methods take a query string for an argument. The format of the string
is described below. In the context of sharing, the query mechanism is beneficial
for the shared-manager process. It is able to perform the query where the data
resides versus the client-process grep locally involving lots of IPC.

 o Basic demonstration

   @keys = $ca->keys( "query string given here" );
   @keys = $ca->keys( "val =~ /pattern/" );

 o Supported operators: =~ !~ eq ne lt le gt ge == != < <= > >=
 o Multiple expressions delimited by :AND or :OR, mixed case allowed

   "key eq 'some key' :or (val > 5 :and val < 9)"
   "key eq some key :or (val > 5 :and val < 9)"
   "key =~ /pattern/i :And val =~ /pattern/i"
   "val eq foo baz :OR key !~ /pattern/i"

   * key matches on keys in the cache
   * likewise, val matches on values

 o Quoting is optional inside the string

   "key =~ /pattern/i :AND val eq 'foo bar'"   # val eq "foo bar"
   "key =~ /pattern/i :AND val eq foo bar"     # val eq "foo bar"

Examples.

 # search capability key/val: =~ !~ eq ne lt le gt ge == != < <= > >=
 # key/val means to match against actual key/val respectively

 @keys  = $ca->keys( "key eq 'some key' :or (val > 5 :and val < 9)" );
 @keys  = $ca->keys( "key eq some key :or (val > 5 :and val < 9)" );

 @keys  = $ca->keys( "key =~ /$pattern/i" );
 @keys  = $ca->keys( "key !~ /$pattern/i" );
 @keys  = $ca->keys( "val =~ /$pattern/i" );
 @keys  = $ca->keys( "val !~ /$pattern/i" );

 %pairs = $ca->pairs( "key == $number" );
 %pairs = $ca->pairs( "key != $number :and val > 100" );
 %pairs = $ca->pairs( "key <  $number :or key > $number" );
 %pairs = $ca->pairs( "val <= $number" );
 %pairs = $ca->pairs( "val >  $number" );
 %pairs = $ca->pairs( "val >= $number" );

 @vals  = $ca->vals( "key eq $string" );
 @vals  = $ca->vals( "key ne $string with space" );
 @vals  = $ca->vals( "key lt $string :or val =~ /$pat1|$pat2/" );
 @vals  = $ca->vals( "val le $string :and val eq 'foo bar'" );
 @vals  = $ca->vals( "val le $string :and val eq foo bar" );
 @vals  = $ca->vals( "val gt $string" );
 @vals  = $ca->vals( "val ge $string" );

=head1 API DOCUMENTATION

This module involves TIE when accessing the object via hash-like behavior.
Both non-shared and shared instances are impacted if doing so. Although likely
fast enough for many use cases, the OO interface is recommended for best
performance.

Accessing an item is likely to involve moving its key to the top of the cache.
Various methods described below state with C<Reorder: Yes> or C<Reorder: No>
as an indication.

The methods C<keys>, C<pairs>, and C<values> return the most frequently
accessed items from the upper section of the cache first before the lower
section. Returned values may not be ordered as expected. This abnormally is
normal for this hybrid LRU-plain implementation. It comes from fetches not
involving LRU movement on keys residing in the upper section of the cache.

When C<max_age> is set, accessing an item which has expired will behave
similarly to a non-existing item.

=over 3

=item new ( { options }, key, value [, key, value, ... ] )

Constructs a new object.

 # non-shared or local construction for use by a single process

 use MCE::Shared::Cache;

 $ca = MCE::Shared::Cache->new(); # max_keys => undef, max_age => undef
 $ca = MCE::Shared::Cache->new( { max_keys => 500 }, @pairs );

 $ca = MCE::Shared::Cache->new( max_keys => "unlimited", max_age => "never" );
 $ca = MCE::Shared::Cache->new( max_keys => undef, max_age => undef ); # ditto
 $ca = MCE::Shared::Cache->new( max_keys => 500, max_age => "1 hour" );
 $ca = MCE::Shared::Cache->new( max_keys => "4 KiB" ); # 4*1024
 $ca = MCE::Shared::Cache->new( max_keys => "1 MiB" ); # 1*1024*1024

 $ca = MCE::Shared::Cache->new( max_age  => "43200 seconds" );
 $ca = MCE::Shared::Cache->new( max_age  => 43200 );   # ditto
 $ca = MCE::Shared::Cache->new( max_age  => "720 minutes" );
 $ca = MCE::Shared::Cache->new( max_age  => "12 hours" );
 $ca = MCE::Shared::Cache->new( max_age  => "0.5 days" );
 $ca = MCE::Shared::Cache->new( max_age  => "1 week" );
 $ca = MCE::Shared::Cache->new( max_age  => "never" ); # or undef
 $ca = MCE::Shared::Cache->new( max_age  => "now" );   # or 0

 $ca->assign( @pairs );

 # construction for sharing with other threads and processes

 use MCE::Shared;

 $ca = MCE::Shared->cache(); # max_keys => undef, max_age => undef
 $ca = MCE::Shared->cache( { max_keys => 500 }, @pairs );

 $ca = MCE::Shared->cache( max_keys => "unlimited", max_age => "never" );
 $ca = MCE::Shared->cache( max_keys => undef, max_age => undef ); # ditto
 $ca = MCE::Shared->cache( max_keys => 500, max_age => "1 hour" );
 $ca = MCE::Shared->cache( max_keys => "4 KiB" ); # 4*1024
 $ca = MCE::Shared->cache( max_keys => "1 MiB" ); # 1*1024*1024

 $ca = MCE::Shared->cache( max_age  => "43200 seconds" );
 $ca = MCE::Shared->cache( max_age  => 43200 );   # ditto
 $ca = MCE::Shared->cache( max_age  => "720 minutes" );
 $ca = MCE::Shared->cache( max_age  => "12 hours" );
 $ca = MCE::Shared->cache( max_age  => "0.5 days" );
 $ca = MCE::Shared->cache( max_age  => "1 week" );
 $ca = MCE::Shared->cache( max_age  => "never" ); # or undef
 $ca = MCE::Shared->cache( max_age  => "now" );   # or 0

 $ca->assign( @pairs );

Reorder: Yes, when given key-value pairs contain duplicate keys

=item assign ( key, value [, key, value, ... ] )

Clears the cache, then sets multiple key-value pairs and returns the number of
keys stored in the cache. This is equivalent to C<clear>, C<mset>.

 $len = $ca->assign( "key1" => "val1", "key2" => "val2" );

Reorder: Yes, when given key-value pairs contain duplicate keys

=item clear

Removes all key-value pairs from the cache.

 $ca->clear;
 %{$ca} = ();

=item delete ( key )

Deletes and returns the value by given key or C<undef> if the key does not
exists in the cache.

 $val = $ca->delete( "some_key" );
 $val = delete $ca->{ "some_key" };

=item del

C<del> is an alias for C<delete>.

=item exists ( key )

Determines if a key exists in the cache.

 if ( $ca->exists( "some_key" ) ) { ... }
 if ( exists $ca->{ "some_key" } ) { ... }

Reorder: No

=item get ( key )

Gets the value of a cache key or C<undef> if the key does not exists.
LRU reordering occurs only if the key is found in the lower section of the
cache. See C<peek> to not promote the key internally to the top of the list.

 $val = $ca->get( "some_key" );
 $val = $ca->{ "some_key" };

Reorder: Yes

=item iterator ( key [, key, ... ] )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns a code reference for iterating a list of key-value pairs stored in
the cache when no arguments are given. Otherwise, returns a code reference for
iterating the given keys in the same order. Keys that do not exist will have
the C<undef> value.

The list of keys to return is set when the closure is constructed. Later keys
added to the hash are not included. Subsequently, the C<undef> value is
returned for deleted keys.

 $iter = $ca->iterator;
 $iter = $ca->iterator( "key1", "key2" );

 while ( my ( $key, $val ) = $iter->() ) {
    ...
 }

Reorder: No

=item iterator ( "query string" )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns a code reference for iterating a list of key-value pairs that match
the given criteria. It returns an empty list if the search found nothing.
The syntax for the C<query string> is described above.

 $iter = $ca->iterator( "val eq some_value" );
 $iter = $ca->iterator( "key eq some_key :AND val =~ /sun|moon|air|wind/" );
 $iter = $ca->iterator( "val eq sun :OR val eq moon :OR val eq foo" );
 $iter = $ca->iterator( "key =~ /$pattern/" );

 while ( my ( $key, $val ) = $iter->() ) {
    ...
 }

Reorder: No

=item keys ( key [, key, ... ] )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns all keys in the cache by most frequently accessed when no arguments
are given. Otherwise, returns the given keys in the same order. Keys that do
not exist will have the C<undef> value. In scalar context, returns the size
of the cache.

 @keys = $ca->keys;
 @keys = $ca->keys( "key1", "key2" );
 $len  = $ca->keys;

Reorder: No

=item keys ( "query string" )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns only keys that match the given criteria. It returns an empty list
if the search found nothing. The syntax for the C<query string> is described
above. In scalar context, returns the size of the resulting list.

 @keys = $ca->keys( "val eq some_value" );
 @keys = $ca->keys( "key eq some_key :AND val =~ /sun|moon|air|wind/" );
 @keys = $ca->keys( "val eq sun :OR val eq moon :OR val eq foo" );
 $len  = $ca->keys( "key =~ /$pattern/" );

Reorder: No

=item len ( key )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns the size of the cache when no arguments are given. For the given key,
returns the length of the value stored at key or the C<undef> value if the
key does not exists.

 $size = $ca->len;
 $len  = $ca->len( "key1" );
 $len  = length $ca->{ "key1" };

Reorder: Yes, only when key is given

=item max_age ( [ secs ] )

Returns the maximum age set on the cache or "never" if not defined internally.
When seconds is given, it adjusts any keys by subtracting or adding the time
difference accordingly.

 $age = $ca->max_age;

 $ca->max_age( "43200 seconds" );
 $ca->max_age( 43200 );     # ditto
 $ca->max_age( "720 minutes" );
 $ca->max_age( "12 hours" );
 $ca->max_age( "0.5 days" );
 $ca->max_age( "1 week" );
 $ca->max_age( "never" );   # or undef
 $ca->max_age( "now" );     # or 0

=item max_keys ( [ size ] )

Returns the size limit set on the cache or "unlimited" if not defined
internally. When size is given, it adjusts the cache accordingly to the
new size by pruning the head of the list if necessary.

 $size = $ca->max_size;

 $ca->max_keys( "unlimited" );
 $ca->max_keys( undef );    # ditto
 $ca->max_keys( "4 KiB" );  # 4*1024
 $ca->max_keys( "1 MiB" );  # 1*1024*1024
 $ca->max_keys( 500 );

=item mdel ( key [, key, ... ] )

Deletes one or more keys in the cache and returns the number of keys deleted.
A given key which does not exist in the cache is not counted.

 $cnt = $ca->mdel( "key1", "key2" );

=item mexists ( key [, key, ... ] )

Returns a true value if all given keys exists in the cache. A false value is
returned otherwise.

 if ( $ca->mexists( "key1", "key2" ) ) { ... }

Reorder: No

=item mget ( key [, key, ... ] )

Gets the values of all given keys. It returns C<undef> for keys which do not
exists in the cache.

 ( $val1, $val2 ) = $ca->mget( "key1", "key2" );

Reorder: Yes

=item mset ( key, value [, key, value, ... ] )

Sets multiple key-value pairs in a cache and returns the number of keys stored
in the cache.

 $len = $ca->mset( "key1" => "val1", "key2" => "val2" );

Reorder: Yes

=item merge

C<merge> is an alias for C<mset>.

=item pairs ( key [, key, ... ] )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns key-value pairs in the cache by most frequently accessed when no
arguments are given. Otherwise, returns key-value pairs for the given keys
in the same order. Keys that do not exist will have the C<undef> value.
In scalar context, returns the size of the cache.

 @pairs = $ca->pairs;
 @pairs = $ca->pairs( "key1", "key2" );
 $len   = $ca->pairs;

Reorder: No

=item pairs ( "query string" )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns only key-value pairs that match the given criteria. It returns an
empty list if the search found nothing. The syntax for the C<query string> is
described above. In scalar context, returns the size of the resulting list.

 @pairs = $ca->pairs( "val eq some_value" );
 @pairs = $ca->pairs( "key eq some_key :AND val =~ /sun|moon|air|wind/" );
 @pairs = $ca->pairs( "val eq sun :OR val eq moon :OR val eq foo" );
 $len   = $ca->pairs( "key =~ /$pattern/" );

Reorder: No

=item peek ( key )

Same as C<get> without changing the order of the keys. Gets the value of a
cache key or C<undef> if the key does not exists.

 $val = $ca->get( "some_key" );
 $val = $ca->{ "some_key" };

Reorder: No

=item pipeline ( [ func1, @args ], [ func2, @args ], ... )

Combines multiple commands for the object to be processed serially. For shared
objects, the call is made atomically due to single IPC to the shared-manager
process. The C<pipeline> method is fully C<wantarray>-aware and receives a list
of commands and their arguments. In scalar or list context, it returns data
from the last command in the pipeline.

 @vals = $ca->pipeline(                     # ( "a_a", "b_b", "c_c" )
    [ "set", foo => "a_a" ],
    [ "set", bar => "b_b" ],
    [ "set", baz => "c_c" ],
    [ "mget", qw/ foo bar baz / ]
 );

 $len = $ca->pipeline(                      # 3, same as $ca->len
    [ "set", foo => "i_i" ],
    [ "set", bar => "j_j" ],
    [ "set", baz => "k_k" ],
    [ "len" ]
 );

 $ca->pipeline(
    [ "set", foo => "m_m" ],
    [ "set", bar => "n_n" ],
    [ "set", baz => "o_o" ]
 );

Reorder: Very likely, see API on given method

=item pipeline_ex ( [ func1, @args ], [ func2, @args ], ... )

Same as C<pipeline>, but returns data for every command in the pipeline.

 @vals = $ca->pipeline_ex(                  # ( "a_a", "b_b", "c_c" )
    [ "set", foo => "a_a" ],
    [ "set", bar => "b_b" ],
    [ "set", baz => "c_c" ]
 );

Reorder: Very likely, see API on given command

=item purge ( )

A utility method for purging any *tombstones* in the keys array. It also
resets a couple counters internally. Expired items are also purged when
max_age is defined.

 $ca->purge;

=item set ( key, value )

Sets the value of the given cache key and returns its new value.

 $val = $ca->set( "key", "value" );
 $val = $ca->{ "key" } = "value";

Reorder: Yes

=item values ( key [, key, ... ] )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns all values in the cache by most frequently accessed when no arguments
are given. Otherwise, returns values for the given keys in the same order.
Keys that do not exist will have the C<undef> value. In scalar context,
returns the size of the cache.

 @vals = $ca->values;
 @vals = $ca->values( "key1", "key2" );
 $len  = $ca->values;

Reorder: No

=item values ( "query string" )

When C<max_age> is set, prunes any expired keys at the head of the list.

Returns only values that match the given criteria. It returns an empty list
if the search found nothing. The syntax for the C<query string> is described
above. In scalar context, returns the size of the resulting list.

 @vals = $ca->values( "val eq some_value" );
 @vals = $ca->values( "key eq some_key :AND val =~ /sun|moon|air|wind/" );
 @vals = $ca->values( "val eq sun :OR val eq moon :OR val eq foo" );
 $len  = $ca->values( "key =~ /$pattern/" );

Reorder: No

=item vals

C<vals> is an alias for C<values>.

=back

=head1 SUGAR METHODS

This module is equipped with sugar methods to not have to call C<set>
and C<get> explicitly. In shared context, the benefit is atomicity and
reduction in inter-process communication.

The API resembles a subset of the Redis primitives
L<http://redis.io/commands#strings> with key representing the cache key.

=over 3

=item append ( key, string )

Appends a value to a key and returns its new length.

 $len = $ca->append( $key, "foo" );

Reorder: Yes

=item decr ( key )

Decrements the value of a key by one and returns its new value.

 $num = $ca->decr( $key );

Reorder: Yes

=item decrby ( key, number )

Decrements the value of a key by the given number and returns its new value.

 $num = $ca->decrby( $key, 2 );

Reorder: Yes

=item getdecr ( key )

Decrements the value of a key by one and returns its old value.

 $old = $ca->getdecr( $key );

Reorder: Yes

=item getincr ( key )

Increments the value of a key by one and returns its old value.

 $old = $ca->getincr( $key );

Reorder: Yes

=item getset ( key, value )

Sets the value of a key and returns its old value.

 $old = $ca->getset( $key, "baz" );

Reorder: Yes

=item incr ( key )

Increments the value of a key by one and returns its new value.

 $num = $ca->incr( $key );

Reorder: Yes

=item incrby ( key, number )

Increments the value of a key by the given number and returns its new value.

 $num = $ca->incrby( $key, 2 );

Reorder: Yes

=back

=head1 PERFORMANCE TESTING

One might want to benchmark this module. If yes, remember to use the non-shared
construction for running on a single core.

 use MCE::Shared::Cache;

 my $cache = MCE::Shared::Cache->new( max_keys => 500_000 );

Otherwise, the following is a parallel version for a L<benchmark script|https://blog.celogeek.com/201401/426/perl-benchmark-cache-with-expires-and-max-size> found on the web. The serial version was created by Celogeek for benchmarking various caching modules.

The MCE C<progress> option makes it possible to track progress while running
among many workers. Being a parallel script means that it involves IPC to and
from the shared-manager process, where the data resides. In regards to IPC,
fetches may take longer on Linux versus running on Darwin or FreeBSD.

 #!/usr/bin/perl

 use strict;
 use warnings;
 use feature 'say', 'state';

 use Digest::MD5 qw/md5_base64/;
 use Time::HiRes qw/time/;
 use Proc::ProcessTable;

 use MCE 1.814;   # requires 1.814 minimally
 use MCE::Shared;

 sub get_current_process_memory {
     state $pt = Proc::ProcessTable->new;
     my %info = map { $_->pid => $_ } @{$pt->table};
     return $info{$$}->rss;
 }

 $| = 1; srand(0);

 # construct shared variables
 # MCE::Shared handles serialization automatically, when needed

 my $c     = MCE::Shared->cache( max_keys => 500_000 );
 my $found = MCE::Shared->scalar( 0 );

 # construct and spawn MCE workers
 # workers increment a local variable $f

 my $mce = MCE->new(
     chunk_size  => 4000,
     max_workers => 4,
     user_func   => sub {
         my ($mce, $chunk_ref, $chunk_id) = @_;
         if ( $mce->user_args()->[0] eq 'setter' ) {
             for ( @{ $chunk_ref } ) { $c->set($_, {md5 => $_})  }
         }
         else {
             my $f = 0;
             for ( @{ $chunk_ref } ) { $f++ if ref $c->get($_) eq 'HASH' }
             $found->incrby($f);
         }
     }
 )->spawn();

 say "Mapping";
 my @todo = map { md5_base64($_) } (1..600_000);
 say "Starting";

 my $mem = get_current_process_memory();
 my ($read, $write);

 {
     my $s = time;
     $mce->process({
         progress  => sub { print "Write: $_[0]\r" },
         user_args => [ 'setter' ],
     }, \@todo);
     $write = time - $s;
 }

 say "Write: ", sprintf("%0.3f", scalar(@todo) / $write);

 {
     my $s = time;
     $found->set(0);
     $mce->process({
         progress  => sub { print "Read $_[0]\r" },
         user_args => [ 'getter' ],
     }, \@todo);
     $read = time - $s;
 }

 $mce->shutdown();

 say "Read : ", sprintf("%0.3f", scalar(@todo) / $read);
 say "Found: ", $found->get();
 say "Mem  : ", get_current_process_memory() - $mem;

The C<progress> option is further described on Metacpan. Several examples
are provided, accommodating all input data-types in MCE.

L<Progress Demonstrations|https://metacpan.org/pod/MCE::Core#PROGRESS-DEMONSTRATIONS>

=head1 SEE ALSO

=over 3

=item * L<CHI>

=item * L<Cache::FastMmap>

=item * L<Cache::LRU>

=item * L<Cache::Ref>

=item * L<Tie::Cache::LRU>

=item * L<Tie::Cache::LRU::Expires>

=back

=head1 INDEX

L<MCE|MCE>, L<MCE::Hobo>, L<MCE::Shared>

=head1 AUTHOR

Mario E. Roy, S<E<lt>marioeroy AT gmail DOT comE<gt>>

=cut

