###############################################################################
## ----------------------------------------------------------------------------
## A pure-Perl in-memory data store.
##
###############################################################################

package MCE::Shared::Minidb;

use strict;
use warnings;

use 5.010001;

no warnings qw( threads recursion uninitialized numeric );

our $VERSION = '1.885';

use MCE::Shared::Base ();
use base 'MCE::Shared::Base::Common';

use MCE::Shared::Ordhash ();
use MCE::Shared::Array ();
use MCE::Shared::Hash ();

use overload (
   q("")    => \&MCE::Shared::Base::_stringify,
   q(0+)    => \&MCE::Shared::Base::_numify,
   fallback => 1
);

sub new {
   # Parallel Hashes: [ HoH, HoA ]
   bless [
      MCE::Shared::Ordhash->new(),  # Hash of Hashes (HoH)
      MCE::Shared::Ordhash->new(),  # Hash of Lists (HoA)
   ], shift;
}

###############################################################################
## ----------------------------------------------------------------------------
## Private methods.
##
###############################################################################

# _hfind ( { getkeys => 1 }, "query string" )
# _hfind ( { getvals => 1 }, "query string" )
# _hfind ( "query string" ) # pairs

sub _hfind {
   my $self   = shift;
   my $params = ref($_[0]) eq 'HASH' ? shift : {};

   if ( @_ == 2 ) {
      my $key = shift;
      return () unless exists($self->[0][0]{ $key });
      $self->[0][0]{ $key }->_find($params, @_);
   }
   else {
      my $query = shift;
      $params->{'hfind'} = undef;

      MCE::Shared::Base::_find_hash(
         $self->[0][0], $params, $query, $self->[0]
      );
   }
}

# _lfind ( { getkeys => 1 }, "query string" )
# _lfind ( { getvals => 1 }, "query string" )
# _lfind ( "query string" ) # pairs

sub _lfind {
   my $self   = shift;
   my $params = ref($_[0]) eq 'HASH' ? shift : {};

   if ( @_ == 2 ) {
      my $key = shift;
      return () unless exists($self->[1][0]{ $key });
      $self->[1][0]{ $key }->_find($params, @_);
   }
   else {
      my $query = shift;
      $params->{'lfind'} = undef;

      MCE::Shared::Base::_find_hash(
         $self->[1][0], $params, $query, $self->[1]
      );
   }
}

# _new_hash ( ) applies to HoH

sub _new_hash {
   MCE::Shared::Hash->new();
}

# _new_list ( ) applies to HoA

sub _new_list {
   MCE::Shared::Array->new();
}

# _qparse ( "select string" )
#
# The select_aref and select_href methods take a select string supporting
# field names or list indices and optionally sort modifiers. The syntax for
# the query string, between :WHERE and :ORDER BY, is the same as described
# in the documentation under the section labeled SYNTAX for QUERY STRING.
#
# The modifiers :WHERE, :AND, :OR, ORDER BY, ASC, DESC, ALPHA may be written
# using mixed case. e.g. :Where
#
# o Hash of Hashes (HoH)
#   "f1 f2 f3 :WHERE f4 > 20 :AND key =~ /foo/ :ORDER BY f5 DESC ALPHA"
#   "f5 f1 f2 :WHERE fN > 40 :AND key =~ /bar/ :ORDER BY key ALPHA"
#   "f5 f1 f2 :WHERE fN > 40 :AND key =~ /bar/"
#   "f5 f1 f2"
#
#    * key matches on keys stored in the primary level hash (H)oH
#
# o Hash of Lists (HoA)
#   "17 15 11 :where 12 > 20 :and key =~ /foo/ :order by 10 desc alpha"
#   "17 15 11 :where 12 > 40 :and key =~ /bar/ :order by key alpha"
#   "17 15 11 :where 12 > 40 :and key =~ /bar/"
#   "17 15 11"
#
#   * key matches on keys stored in the primary level hash (H)oA
#   * above, list indices are given as 17, 15, 11, 12, and 10
#   * the shorter form is allowed e.g. "4 > 20 :AND key =~ /baz/"

sub _qparse {
   my ( $q ) = @_;
   my ( $f, $w, $o );

   if ( $q =~ /^([\S ]*):where[ ]+(.+):order by[ ]+(.+)/i ) {
      ( $f, $w, $o ) = ( $1, $2, $3 );
   }
   elsif ( $q =~ /^([\S ]*):where[ ]+(.+)/i ) {
      ( $f, $w ) = ( $1, $2 );
   }
   elsif ( $q =~ /^([\S ]*):order by[ ]+(.+)/i ) {
      ( $f, $o ) = ( $1, $2 );
   }
   elsif ( $q =~ /^((?:key|\S+)[ ]+(?:=|!|<|>|e|n|l|g)\S?[ ]+\S.*)/ ) {
      ( $w ) = ( $1 );
   }
   elsif ( $q =~ /^([\S ]*)/ ) {
      ( $f ) = ( $1 );
   }

   $f =~ s/[ ]+$//, $w =~ s/[ ]+$//, $o =~ s/[ ]+$//;

   return ( $f, $w, $o );
}

# _hselect_aref ( "select string" ), see _qparse for description
#
#  returns an array containing [ key, aref ] pairs

sub _hselect_aref {
   my ( $self, $query ) = @_;
   my ( $f, $w, $o ) = _qparse($query);

   my @fields = split(' ', $f);
   my $data   = $self->[0][0];

   unless ( @fields ) {
      warn("_hselect_aref: must specify fieldname(s)");
      return ();
   }

   if ( length $w ) {
      my %match = map { $_ => 1 } ( $self->hkeys($w) );
      map { !exists $match{$_} ? () : do {
               my ( $k, @ret ) = ( $_ );
               push @ret, $data->{$k}{$_} for @fields;
               [ $k, \@ret ];
            };
          } ( length $o ? $self->hsort($o) : $self->hkeys() );
   }
   else {
      map { my ( $k, @ret ) = ( $_ );
            push @ret, $data->{$k}{$_} for @fields;
            [ $k, \@ret ];
          } ( length $o ? $self->hsort($o) : $self->hkeys() );
   }
}

# _hselect_href ( "select string" ), see _qparse for description
#
#  returns an array containing [ key, href ] pairs

sub _hselect_href {
   my ( $self, $query ) = @_;
   my ( $f, $w, $o ) = _qparse($query);

   my @fields = split(' ', $f);
   my $data   = $self->[0][0];

   if ( length $w ) {
      my %match = map { $_ => 1 } ( $self->hkeys($w) );
      if ( @fields ) {
         map { !exists $match{$_} ? () : do {
                  my ( $k, %ret ) = ( $_ );
                  $ret{$_} = $data->{$k}{$_} for @fields;
                  [ $k, \%ret ];
               };
             } ( length $o ? $self->hsort($o) : $self->hkeys() );
      }
      else {
         map { !exists $match{$_} ? () : [ $_, { %{ $data->{$_} } } ];
             } ( length $o ? $self->hsort($o) : $self->hkeys() );
      }
   }
   else {
      if ( @fields ) {
         map { my ( $k, %ret ) = ( $_ );
               $ret{$_} = $data->{$k}{$_} for @fields;
               [ $k, \%ret ];
             } ( length $o ? $self->hsort($o) : $self->hkeys() );
      }
      else {
         map { [ $_, { %{ $data->{$_} } } ];
             } ( length $o ? $self->hsort($o) : $self->hkeys() );
      }
   }
}

# _lselect_aref ( "select string" ), see _qparse for description
#
#  returns an array containing [ key, aref ] pairs

sub _lselect_aref {
   my ( $self, $query ) = @_;
   my ( $f, $w, $o ) = _qparse($query);

   my @fields = split(' ', $f);
   my $data   = $self->[1][0];

   if ( length $w ) {
      my %match = map { $_ => 1 } ( $self->lkeys($w) );
      if ( @fields ) {
         map { !exists $match{$_} ? () : do {
                  my ( $k, @ret ) = ( $_ );
                  push @ret, $data->{$k}[$_] for @fields;
                  [ $k, \@ret ];
               };
             } ( length $o ? $self->lsort($o) : $self->lkeys() );
      }
      else {
         map { !exists $match{$_} ? () : [ $_, [ @{ $data->{$_} } ] ];
             } ( length $o ? $self->lsort($o) : $self->lkeys() );
      }
   }
   else {
      if ( @fields ) {
         map { my ( $k, @ret ) = ( $_ );
               push @ret, $data->{$k}[$_] for @fields;
               [ $k, \@ret ];
             } ( length $o ? $self->lsort($o) : $self->lkeys() );
      }
      else {
         map { [ $_, [ @{ $data->{$_} } ] ];
             } ( length $o ? $self->lsort($o) : $self->lkeys() );
      }
   }
}

# _lselect_href ( "select string" ), see _qparse for description
#
#  returns an array containing [ key, href ] pairs

sub _lselect_href {
   my ( $self, $query ) = @_;
   my ( $f, $w, $o ) = _qparse($query);

   my @fields = split(' ', $f);
   my $data = $self->[1][0];

   if ( length $w ) {
      my %match = map { $_ => 1 } ( $self->lkeys($w) );
      if ( @fields ) {
         map { !exists $match{$_} ? () : do {
                  my ( $k, %ret ) = ( $_ );
                  $ret{$_} = $data->{$k}[$_] foreach @fields;
                  [ $k, \%ret ];
               };
             } ( length $o ? $self->lsort($o) : $self->lkeys() );
      }
      else {
         map { !exists $match{$_} ? () : do {
                  my ( $k, %ret ) = ( $_ );
                  $ret{$_} = $data->{$k}[$_] for 0 .. @{ $data->{$k} } - 1;
                  [ $k, \%ret ];
               };
             } ( length $o ? $self->lsort($o) : $self->lkeys() );
      }
   }
   else {
      if ( @fields ) {
         map { my ( $k, %ret ) = ( $_ );
               $ret{$_} = $data->{$k}[$_] foreach @fields;
               [ $k, \%ret ];
             } ( length $o ? $self->lsort($o) : $self->lkeys() );
      }
      else {
         map { my ( $k, %ret ) = ( $_ );
               $ret{$_} = $data->{$k}[$_] for 0 .. @{ $data->{$k} } - 1;
               [ $k, \%ret ];
             } ( length $o ? $self->lsort($o) : $self->lkeys() );
      }
   }
}

# _sort ( HoH, 0, "BY key   [ ASC | DESC ] [ ALPHA ]" )
# _sort ( HoH, 0, "BY field [ ASC | DESC ] [ ALPHA ]" ) e.g. BY address
# _sort ( HoA, 1, "BY key   [ ASC | DESC ] [ ALPHA ]" )
# _sort ( HoA, 1, "BY index [ ASC | DESC ] [ ALPHA ]" ) e.g. BY 9

sub _sort {
   my ( $o, $is_list, $request ) = @_;

   return () unless ( length $request );
   $request =~ s/^[ ]*\bby\b[ ]*//i;

   if ( $request =~ /^[ ]*(\S+)[ ]*(.*)/ ) {
      my ( $f, $modifiers, $alpha, $desc ) = ( $1, $2, 0, 0 );

      $alpha = 1 if ( $modifiers =~ /\balpha\b/i );
      $desc  = 1 if ( $modifiers =~ /\bdesc\b/i );

      # Return sorted keys, leaving the data intact.

      if ( defined wantarray ) {
         if ( $f eq 'key' ) {                         # by key
            if ( $alpha ) { ( $desc )
             ? sort { $b cmp $a } $o->keys
             : sort { $a cmp $b } $o->keys;
            }
            else { ( $desc )
             ? sort { $b <=> $a } $o->keys
             : sort { $a <=> $b } $o->keys;
            }
         }
         else {                                       # by field
            my $d = $o->[0];
            if ( $is_list ) {
               if ( $alpha ) { ( $desc )
                ? sort { $d->{$b}[$f] cmp $d->{$a}[$f] } $o->keys
                : sort { $d->{$a}[$f] cmp $d->{$b}[$f] } $o->keys;
               }
               else { ( $desc )
                ? sort { $d->{$b}[$f] <=> $d->{$a}[$f] } $o->keys
                : sort { $d->{$a}[$f] <=> $d->{$b}[$f] } $o->keys;
               }
            }
            else {
               if ( $alpha ) { ( $desc )
                ? sort { $d->{$b}{$f} cmp $d->{$a}{$f} } $o->keys
                : sort { $d->{$a}{$f} cmp $d->{$b}{$f} } $o->keys;
               }
               else { ( $desc )
                ? sort { $d->{$b}{$f} <=> $d->{$a}{$f} } $o->keys
                : sort { $d->{$a}{$f} <=> $d->{$b}{$f} } $o->keys;
               }
            }
         }
      }

      # Sort in-place otherwise, in void context.

      elsif ( $f eq 'key' ) {                         # by key
         if ( $alpha ) { ( $desc )
          ? $o->_reorder( sort { $b cmp $a } $o->keys )
          : $o->_reorder( sort { $a cmp $b } $o->keys );
         }
         else { ( $desc )
          ? $o->_reorder( sort { $b <=> $a } $o->keys )
          : $o->_reorder( sort { $a <=> $b } $o->keys );
         }
      }
      else {                                          # by field
         my $d = $o->[0];
         if ( $is_list ) {
            if ( $alpha ) { ( $desc )
             ? $o->_reorder( sort { $d->{$b}[$f] cmp $d->{$a}[$f] } $o->keys )
             : $o->_reorder( sort { $d->{$a}[$f] cmp $d->{$b}[$f] } $o->keys );
            }
            else { ( $desc )
             ? $o->_reorder( sort { $d->{$b}[$f] <=> $d->{$a}[$f] } $o->keys )
             : $o->_reorder( sort { $d->{$a}[$f] <=> $d->{$b}[$f] } $o->keys );
            }
         }
         else {
            if ( $alpha ) { ( $desc )
             ? $o->_reorder( sort { $d->{$b}{$f} cmp $d->{$a}{$f} } $o->keys )
             : $o->_reorder( sort { $d->{$a}{$f} cmp $d->{$b}{$f} } $o->keys );
            }
            else { ( $desc )
             ? $o->_reorder( sort { $d->{$b}{$f} <=> $d->{$a}{$f} } $o->keys )
             : $o->_reorder( sort { $d->{$a}{$f} <=> $d->{$b}{$f} } $o->keys );
            }
         }
      }
   }
   else {
      ();
   }
}

###############################################################################
## ----------------------------------------------------------------------------
## Common methods.
##
###############################################################################

# dump ( "file.dat" )

sub dump {
   my ( $self, $file ) = @_;

   if ( length $file ) {
      require Storable unless $INC{'Storable.pm'};

      # purge tombstones
      $self->[0]->purge(), $self->[1]->purge();

      local $@; local $SIG{__DIE__};
      eval { Storable::nstore($self, $file) };

      warn($@), return if $@;
   }
   else {
      warn('Usage: $obj->dump("file.dat")');
      return;
   }

   1;
}

# restore ( "file.dat" )

sub restore {
   my ( $self, $file ) = @_;

   if ( length $file ) {
      require Storable unless $INC{'Storable.pm'};

      local $@; local $SIG{__DIE__};
      my $obj = eval { Storable::retrieve($file) };
      warn($@), return if $@;

      if ( ref($obj) ne 'MCE::Shared::Minidb' ) {
         warn("$file isn't serialized Minidb data: ".ref($obj));
         return;
      }
      $self->[1]->clear(), $self->[1] = delete $obj->[1];
      $self->[0]->clear(), $self->[0] = delete $obj->[0];
   }
   else {
      warn('Usage: $obj->restore("file.dat")');
      return;
   }

   1;
}

# select_aref ( ":lists", "select string" )
# select_aref ( ":hashes", "select string" )
# select_aref ( "select string" )  same as ":hashes"

sub select_aref {
   my ( $self, @query ) = @_;

   if ( $query[0] =~ /^:lists$/i ) {
      shift @query;
      $self->_lselect_aref($query[0]);
   }
   else {
      shift @query if ( $query[0] =~ /^:hashes$/i );
      $self->_hselect_aref($query[0]);
   }
}

# select_href ( ":lists", "select string" )
# select_href ( ":hashes", "select string" )
# select_href ( "select string" )  same as ":hashes"

sub select_href {
   my ( $self, @query ) = @_;

   if ( $query[0] =~ /^:lists$/i ) {
      shift @query;
      $self->_lselect_href($query[0]);
   }
   else {
      shift @query if ( $query[0] =~ /^:hashes$/i );
      $self->_hselect_href($query[0]);
   }
}

###############################################################################
## ----------------------------------------------------------------------------
## Hash of Hashes (HoH).
##
###############################################################################

# hset ( key, field, value [, field, value, ... ] )

sub hset {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   if ( @_ ) {
      $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
      if ( @_ == 2 ) {
         $self->[0][0]{ $key }{ $_[0] } = $_[1];
      } else {
         $self->[0][0]{ $key }->mset(@_);
      }
   }
   else {
      return;
   }
}

# hget ( key, field [, field, ... ] )
# hget ( key )

sub hget {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   if ( @_ ) {
      return unless exists($self->[0][0]{ $key });
      if ( @_ == 1 ) {
         $self->[0][0]{ $key }{ $_[0] };
      } else {
         $self->[0][0]{ $key }->mget(@_);
      }
   }
   else {
      $self->[0][0]{ $key };
   }
}

# hdel ( key, field [, field, ... ] )
# hdel ( key )

sub hdel {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   if ( @_ ) {
      return unless exists($self->[0][0]{ $key });
      if ( @_ == 1 ) {
         delete $self->[0][0]{ $key }{ $_[0] };
      } else {
         $self->[0][0]{ $key }->mdel(@_);
      }
   }
   else {
      $self->[0]->del($key);
   }
}

# hexists ( key, field [, field, ... ] )
# hexists ( key )

sub hexists {
   my ( $self, $key ) = ( shift, shift );
   return '' unless length($key);
   if ( @_ ) {
      return '' unless exists($self->[0][0]{ $key });
      if ( @_ == 1 ) {
         exists $self->[0][0]{ $key }{ $_[0] };
      } else {
         $self->[0][0]{ $key }->mexists(@_);
      }
   }
   else {
      exists $self->[0][0]{ $key };
   }
}

# hclear ( key )
# hclear ( )

sub hclear {
   my ( $self, $key ) = @_;
   if ( @_ > 1 ) {
      return unless exists($self->[0][0]{ $key });
      $self->[0][0]{ $key }->clear();
   }
   else {
      $self->[0]->clear();
   }
}

# hkeys ( key, [ field [, field, ... ] ] )
# hkeys ( key, "query string" )
# hkeys ( "query string" )
# hkeys ( )

sub hkeys {
   my $self = shift;

   if ( @_ == 1 && $_[0] =~ /^(?:key|\S+)[ ]+\S\S?[ ]+\S/ ) {
      $self->_hfind({ getkeys => 1 }, @_);
   }
   elsif ( @_ ) {
      my $key = shift;
      return () unless exists($self->[0][0]{ $key });
      $self->[0][0]{ $key }->keys(@_);
   }
   else {
      $self->[0]->keys();
   }
}

# hpairs ( key, [ field [, field, ... ] ] )
# hpairs ( key, "query string" )
# hpairs ( "query string" )
# hpairs ( )

sub hpairs {
   my $self = shift;

   if ( @_ == 1 && $_[0] =~ /^(?:key|\S+)[ ]+\S\S?[ ]+\S/ ) {
      $self->_hfind({}, @_);
   }
   elsif ( @_ ) {
      my $key = shift;
      return () unless exists($self->[0][0]{ $key });
      $self->[0][0]{ $key }->pairs(@_);
   }
   else {
      $self->[0]->pairs();
   }
}

# hvals ( key, [ field [, field, ... ] ] )
# hvals ( key, "query string" )
# hvals ( "query string" )
# hvals ( )

sub hvals {
   my $self = shift;

   if ( @_ == 1 && $_[0] =~ /^(?:key|\S+)[ ]+\S\S?[ ]+\S/ ) {
      $self->_hfind({ getvals => 1 }, @_);
   }
   elsif ( @_ ) {
      my $key = shift;
      return () unless exists($self->[0][0]{ $key });
      $self->[0][0]{ $key }->vals(@_);
   }
   else {
      $self->[0]->vals();
   }
}

# hshift ( )

sub hshift {
   $_[0]->[0]->shift();
}

# hsort ( "BY key   [ ASC | DESC ] [ ALPHA ]" )
# hsort ( "BY field [ ASC | DESC ] [ ALPHA ]" )

sub hsort {
   my ( $self, $request ) = @_;
   return () unless ( @_ == 2 );
   _sort($self->[0], 0, $request);
}

# happend ( key, field, string )

sub happend {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   length( $self->[0][0]{ $key }{ $_[2] } .= $_[3] // '' );
}

# hassign ( key, field, value [, field, value, ... ] )

sub hassign {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   $self->[0][0]{ $key }->assign(@_);
}

# hdecr ( key, field )

sub hdecr {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   --$self->[0][0]{ $key }{ $_[2] };
}

# hdecrby ( key, field, number )

sub hdecrby {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   $self->[0][0]{ $key }{ $_[2] } -= $_[3] || 0;
}

# hincr ( key, field )

sub hincr {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   ++$self->[0][0]{ $key }{ $_[2] };
}

# hincrby ( key, field, number )

sub hincrby {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   $self->[0][0]{ $key }{ $_[2] } += $_[3] || 0;
}

# hgetdecr ( key, field )

sub hgetdecr {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   $self->[0][0]{ $key }{ $_[2] }-- // 0;
}

# hgetincr ( key, field )

sub hgetincr {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   $self->[0][0]{ $key }{ $_[2] }++ // 0;
}

# hgetset ( key, field, value )

sub hgetset {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   $self->[0][0]{ $key }->getset(@_);
}

# hsetnx ( key, field, value )

sub hsetnx {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   $self->[0]->set($key, _new_hash()) unless exists($self->[0][0]{ $key });
   $self->[0][0]{ $key }->setnx(@_);
}

# hlen ( key, field )
# hlen ( key )
# hlen ( )

sub hlen {
   my $self = shift;
   if ( @_ ) {
      my $key = shift;
      return 0 unless exists($self->[0][0]{ $key });
      $self->[0][0]{ $key }->len(@_);
   }
   else {
      $self->[0]->len();
   }
}

###############################################################################
## ----------------------------------------------------------------------------
## Hash of Lists (HoA).
##
###############################################################################

# lset ( key, index, value [, index, value, ... ] )

sub lset {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   if ( @_ ) {
      $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
      if ( @_ == 2 ) {
         $self->[1][0]{ $key }[ $_[0] ] = $_[1];
      } else {
         $self->[1][0]{ $key }->mset(@_);
      }
   }
   else {
      return;
   }
}

# lget ( key, index [, index, ... ] )
# lget ( key )

sub lget {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   if ( @_ ) {
      return unless exists($self->[1][0]{ $key });
      if ( @_ == 1 ) {
         $self->[1][0]{ $key }[ $_[0] ];
      } else {
         $self->[1][0]{ $key }->mget(@_);
      }
   }
   else {
      $self->[1][0]{ $key };
   }
}

# ldel ( key, index [, index, ... ] )
# ldel ( key )

sub ldel {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   if ( @_ ) {
      return unless exists($self->[1][0]{ $key });
      if ( @_ == 1 ) {
         delete $self->[1][0]{ $key }[ $_[0] ];
      } else {
         $self->[1][0]{ $key }->mdel(@_);
      }
   }
   else {
      $self->[1]->del($key);
   }
}

# lexists ( key, index [, index, ... ] )
# lexists ( key )

sub lexists {
   my ( $self, $key ) = ( shift, shift );
   return '' unless length($key);
   if ( @_ ) {
      return '' unless exists($self->[1][0]{ $key });
      if ( @_ == 1 ) {
         exists $self->[1][0]{ $key }[ $_[0] ];
      } else {
         $self->[1][0]{ $key }->mexists(@_);
      }
   }
   else {
      exists $self->[1][0]{ $key };
   }
}

# lclear ( key )
# lclear ( )

sub lclear {
   my ( $self, $key ) = @_;
   if ( @_ > 1 ) {
      return unless exists($self->[1][0]{ $key });
      $self->[1][0]{ $key }->clear();
   }
   else {
      $self->[1]->clear();
   }
}

# lrange ( key, start, stop )

sub lrange {
   my ( $self, $key ) = ( shift, shift );
   return () unless length($key) && exists($self->[1][0]{ $key });
   $self->[1][0]{ $key }->range(@_);
}

# lsplice ( key, offset [, length [, list ] ] )

sub lsplice {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key) && scalar(@_);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   $self->[1][0]{ $key }->splice(@_);
}

# lpop ( key )

sub lpop {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key) && exists($self->[1][0]{ $key });
   shift @{ $self->[1][0]{ $key } };
}

# lpush ( key, value [, value, ... ] )

sub lpush {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key) && scalar(@_);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   unshift @{ $self->[1][0]{ $key } }, @_;
}

# rpop ( key )

sub rpop {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key) && exists($self->[1][0]{ $key });
   pop @{ $self->[1][0]{ $key } };
}

# rpush ( key, value [, value, ... ] )

sub rpush {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key) && scalar(@_);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   push @{ $self->[1][0]{ $key } }, @_;
}

# lkeys ( key, [ index [, index, ... ] ] )
# lkeys ( key, "query string" )
# lkeys ( "query string" )
# lkeys ( )

sub lkeys {
   my $self = shift;

   if ( @_ == 1 && $_[0] =~ /^(?:key|\S+)[ ]+\S\S?[ ]+\S/ ) {
      $self->_lfind({ getkeys => 1 }, @_);
   }
   elsif ( @_ ) {
      my $key = shift;
      return () unless exists($self->[1][0]{ $key });
      $self->[1][0]{ $key }->keys(@_);
   }
   else {
      $self->[1]->keys();
   }
}

# lpairs ( key, [ index [, index, ... ] ] )
# lpairs ( key, "query string" )
# lpairs ( "query string" )
# lpairs ( )

sub lpairs {
   my $self = shift;

   if ( @_ == 1 && $_[0] =~ /^(?:key|\S+)[ ]+\S\S?[ ]+\S/ ) {
      $self->_lfind({}, @_);
   }
   elsif ( @_ ) {
      my $key = shift;
      return () unless exists($self->[1][0]{ $key });
      $self->[1][0]{ $key }->pairs(@_);
   }
   else {
      $self->[1]->pairs();
   }
}

# lvals ( key, [ index [, index, ... ] ] )
# lvals ( key, "query string" )
# lvals ( "query string" )
# lvals ( )

sub lvals {
   my $self = shift;

   if ( @_ == 1 && $_[0] =~ /^(?:key|\S+)[ ]+\S\S?[ ]+\S/ ) {
      $self->_lfind({ getvals => 1 }, @_);
   }
   elsif ( @_ ) {
      my $key = shift;
      return () unless exists($self->[1][0]{ $key });
      $self->[1][0]{ $key }->vals(@_);
   }
   else {
      $self->[1]->vals();
   }
}

# lshift ( )

sub lshift {
   $_[0]->[1]->shift();
}

# lsort ( "BY key   [ ASC | DESC ] [ ALPHA ]" )
# lsort ( "BY index [ ASC | DESC ] [ ALPHA ]" )
#
# lsort ( key, "BY key [ ASC | DESC ] [ ALPHA ]" )
# lsort ( key, "BY val [ ASC | DESC ] [ ALPHA ]" )

sub lsort {
   my ( $self, $arg1, $arg2 ) = @_;
   if ( @_ == 2 ) {
      _sort($self->[1], 1, $arg1);
   }
   else {
      return () unless ( @_ == 3 && exists($self->[1][0]{ $arg1 }) );
      $self->[1][0]{ $arg1 }->sort($arg2);
   }
}

# lappend ( key, index, string )

sub lappend {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   length( $self->[1][0]{ $key }[ $_[2] ] .= $_[3] // '' );
}

# lassign ( key, value [, value, ... ] )

sub lassign {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   $self->[1][0]{ $key }->assign(@_);
}

# ldecr ( key, index )

sub ldecr {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   --$self->[1][0]{ $key }[ $_[2] ];
}

# ldecrby ( key, index, number )

sub ldecrby {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   $self->[1][0]{ $key }[ $_[2] ] -= $_[3] || 0;
}

# lincr ( key, index )

sub lincr {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   ++$self->[1][0]{ $key }[ $_[2] ];
}

# lincrby ( key, index, number )

sub lincrby {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   $self->[1][0]{ $key }[ $_[2] ] += $_[3] || 0;
}

# lgetdecr ( key, index )

sub lgetdecr {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   $self->[1][0]{ $key }[ $_[2] ]-- // 0;
}

# lgetincr ( key, index )

sub lgetincr {
   my ( $self, $key ) = @_;
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   $self->[1][0]{ $key }[ $_[2] ]++ // 0;
}

# lgetset ( key, index, value )

sub lgetset {
   my ( $self, $key ) = ( shift, shift );
   return unless length($key);
   $self->[1]->set($key, _new_list()) unless exists($self->[1][0]{ $key });
   $self->[1][0]{ $key }->getset(@_);
}

# llen ( key, index )
# llen ( key )
# llen ( )

sub llen {
   my $self = shift;
   if ( @_ ) {
      my $key = shift;
      return 0 unless exists($self->[1][0]{ $key });
      $self->[1][0]{ $key }->len(@_);
   }
   else {
      $self->[1]->len();
   }
}

1;

__END__

###############################################################################
## ----------------------------------------------------------------------------
## Module usage.
##
###############################################################################

=head1 NAME

MCE::Shared::Minidb - A pure-Perl in-memory data store

=head1 VERSION

This document describes MCE::Shared::Minidb version 1.885

=head1 DESCRIPTION

A tiny in-memory NoSQL-like database for use as a standalone or managed by
L<MCE::Shared>.

This module was created mainly for having an efficient manner in which to
manipulate hashes-of-hashes (HoH) and hashes-of-lists (HoA) structures with
MCE::Shared. An application may choose to use both structures or one or the
other. Internally, both structures reside in memory simultaneously.

 sub new {
    # Dual top-level hashes: [ HoH, HoA ]
    bless [
       MCE::Shared::Ordhash->new(),  # Stores Hash-of-Hashes (HoH)
       MCE::Shared::Ordhash->new(),  # Stores Hash-of-Lists  (HoA)
    ], shift;
 }

 # each Ho(H) key => MCE::Shared::Hash->new()
 # each Ho(A) key => MCE::Shared::Array->new()

Several methods described below may resemble the C<Redis> API. It is not the
intent for this module to become 100% compatible.

=head1 SYNOPSIS

 # non-shared or local construction for use by a single process

 use MCE::Shared::Minidb;

 my $db = MCE::Shared::Minidb->new();

 # construction for sharing with other threads and processes

 use MCE::Shared;

 my $db = MCE::Shared->minidb();

 # Hash of Hashes (HoH)
 # Update the hash stored at key1/key2

 $db->hset( "key1", "f1", "foo" );
 $db->hset( "key2", "f1", "bar", "f2", "baz" );

 $val = $db->hget( "key2", "f2" );  # "baz"

 # Hash of Lists (HoA)
 # Update the list stored at key1/key2

 $db->lset( "key1", 0, "foo" );
 $db->lset( "key2", 0, "bar", 1, "baz" );

 $val = $db->lget( "key2", 1 );     # "baz"

 # pipeline, provides atomicity for shared objects, MCE::Shared v1.09+

 $db->pipeline(
    [ "lset", "key1", 0, "foo" ],
    [ "hset", "key1", "field1", "foo" ],
    [ "lset", "key2", 0, "bar", 1, "baz" ],
    [ "hset", "key2", "field1", "bar", "field2", "baz" ]
 );

=head1 SYNTAX for QUERY STRING

Several methods take a query string for an argument. The format of the string
is described below. In the context of sharing, the query mechanism is beneficial
for the shared-manager process. It is able to perform the query where the data
resides versus the client-process grep locally involving lots of IPC.

 o Basic demonstration

   # query the hash stored at "some key"
   @keys = $db->hkeys( "some key", "query string given here" );
   @keys = $db->hkeys( "some key", "val =~ /pattern/" );

   # query the top-level hash (H)oH
   @keys = $db->hkeys( "query string given here" );
   @keys = $db->hkeys( "key =~ /pattern/" );

 o Supported operators: =~ !~ eq ne lt le gt ge == != < <= > >=
 o Multiple expressions delimited by :AND or :OR, mixed case allowed

   "key eq 'some key' :or (field > 5 :and field < 9)"
   "key eq some key :or (field > 5 :and field < 9)"
   "key =~ /pattern/i :And field =~ /pattern/i"   # HoH
   "key =~ /pattern/i :And index =~ /pattern/i"   # HoA
   "index eq foo baz :OR key !~ /pattern/i"       # e.g. 9 eq "foo baz"

   * key   matches on primary keys in the hash (H)oH or (H)oA
   * field matches on HoH->{key}{field} e.g. address
   * index matches on HoA->{key}[index] e.g. 9

 o Quoting is optional inside the string
 o Primary keys (H)oH may have spaces, but not so for field_names

   "key =~ /pattern/i :AND field eq 'foo bar'"    # address eq "foo bar"
   "key =~ /pattern/i :AND field eq foo bar"      # address eq "foo bar"

   "key =~ 'some key' :AND 'some_field' eq 'foo bar'"  # ok: some_field
   "key =~ some key :AND some_field eq foo bar"

   "key =~ 'some key' :AND 'some field' eq 'foo bar'"  # fail: some field
   "key =~ some key :AND some field eq foo bar"

Examples.

 # search capability key/val: =~ !~ eq ne lt le gt ge == != < <= > >=
 # key/val means to match against actual key/val respectively

 # a query made to a hash stored at "some key"
 # note that "fieldNames" must not have spaces

   @keys  = $db->hkeys(
      "some key", "key eq some_field :or (val > 5 :and val < 9)"
   );

 # queries made to a list stored at "some key"

   @pairs = $db->lpairs( "some key", "key >= 50 :AND val =~ /sun|moon/" );
   @pairs = $db->lpairs( "some key", "val eq sun :OR val eq moon" );

 # the key modifier is the only thing possible for top-level queries
 # reason: value equals a hash or list reference containing 2nd-level data

   @keys  = $db->hkeys( "key eq 'some key'" );
   @keys  = $db->hkeys( "key eq some key" );

   @keys  = $db->hkeys( "key =~ /$pattern/i" );
   @keys  = $db->hkeys( "key !~ /$pattern/i" );

   %pairs = $db->hpairs( "key == $number" );
   %pairs = $db->hpairs( "key != $number" );
   %pairs = $db->hpairs( "key <  $number :or key > $number" );
   %pairs = $db->hpairs( "key <= $number" );
   %pairs = $db->hpairs( "key >  $number" );
   %pairs = $db->hpairs( "key >= $number" );

   @vals  = $db->hvals( "key eq $string" );
   @vals  = $db->hvals( "key ne $string with space" );
   @vals  = $db->hvals( "key lt $string :or key =~ /$pat1|$pat2/" );
   @vals  = $db->hvals( "key le $string :or key eq 'foo bar'" );
   @vals  = $db->hvals( "key le $string :or key eq foo bar" );
   @vals  = $db->hvals( "key gt $string" );
   @vals  = $db->hvals( "key ge $string" );

 # see select_aref and select_href below for db-like queries against
 # the underlying Ho(A) and Ho(H) structures respectively

=head1 API DOCUMENTATION - DB

=head2 MCE::Shared::Minidb->new ()

=head2 MCE::Shared->minidb ()

Constructs an empty in-memory C<HoH> and C<HoA> key-store database structure.

 # non-shared or local construction for use by a single process

 use MCE::Shared::Minidb;

 $db = MCE::Shared::Minidb->new();

 # construction for sharing with other threads and processes

 use MCE::Shared;

 $db = MCE::Shared->minidb();

=head2 dump ( "file.dat" )

Dumps the in-memory content to a file.

 $db->dump( "content.dat" );

=head2 pipeline ( [ func1, @args ], [ func2, @args ], ... )

Combines multiple commands for the object to be processed serially. For shared
objects, the call is made atomically due to single IPC to the shared-manager
process. The C<pipeline> method is fully C<wantarray>-aware and receives a list
of commands and their arguments. In scalar or list context, it returns data
from the last command in the pipeline.

 @vals = $db->pipeline(                     # ( "bar", "baz" )
    [ "hset", "some_key", "f1", "bar", "f2", "baz" ],
    [ "hget", "some_key", "f1", "f2" ]
 );

 $len = $db->pipeline(                      # 2, same as $db->hlen("key2)
    [ "hset", "some_key", "f1", "bar", "f2", "baz" ],
    [ "hlen", "some_key" ]
 );

 $db->pipeline(
    [ "lset", "key1", 0, "foo" ],
    [ "hset", "key1", "field1", "foo" ],
    [ "lset", "key2", 0, "bar", 1, "baz" ],
    [ "hset", "key2", "field1", "bar", "field2", "baz" ]
 );

Current API available since 1.809.

=head2 pipeline_ex ( [ func1, @args ], [ func2, @args ], ... )

Same as C<pipeline>, but returns data for every command in the pipeline.

 @vals = $db->pipeline_ex(                  # ( "bar", "baz" )
    [ "hset", "key3", "field1", "bar" ],
    [ "hset", "key3", "field2", "baz" ]
 );

 $chunk_size = 3;

 $db->hset("some_key", chunk_id => 0);
 $db->lassign("some_key", @ARGV);

 while (1) {
    ($chunk_id, @next) = $db->pipeline_ex(
       [ "hincr",   "some_key", "chunk_id"     ],
       [ "lsplice", "some_key", 0, $chunk_size ]
    );

    last unless @next;

    ...
 }

Current API available since 1.809.

=head2 restore ( "file.dat" )

Restores the in-memory content from a file.

 $db->restore( "content.dat" );

=head2 select_aref ( ":hashes", "select string" )

=head2 select_href ( ":hashes", "select string" )

Returns a list containing C<[ key, aref ]> pairs or C<[ key, href ]> pairs
from the hash of hashes (HoH).

The C<select_aref> and C<select_href> methods take a select string supporting
field names and optionally sort modifiers. The syntax for the query string,
between C<:WHERE> and C<:ORDER BY>, is the same as described above.

The modifiers C<:WHERE>, C<:AND>, C<:OR>, C<ORDER BY>, C<ASC>, C<DESC>, and
C<ALPHA> may be mixed case. e.g. C<:Where>

HoH select string:

 "f1 f2 f3 :WHERE f4 > 20 :AND key =~ /foo/ :ORDER BY f5 DESC ALPHA"
 "f5 f1 f2 :WHERE fN > 40 :AND key =~ /bar/ :ORDER BY key ALPHA"
 "f5 f1 f2 :WHERE fN > 40 :AND key =~ /bar/"
 "f5 f1 f2"

 * key matches on keys stored in the primary level hash (H)oH

HoH select synopsis:

 @rows = $db->select_aref( ":hashes", "some_field :ORDER BY some_field" );
 @rows = $db->select_aref( ":hashes", "f2 f6 f5 :ORDER BY f4" );

 # $rows[0] = [ "key1", [ "f2_val", "f6_val", "f5_val" ] ]
 # $rows[1] = [ "key2", [ "f2_val", "f6_val", "f5_val" ] ]
 # $rows[2] = [ "key3", [ "f2_val", "f6_val", "f5_val" ] ]
 # ...
 # $rows[N] = [ "keyN", [ "f2_val", "f6_val", "f5_val" ] ]

 @rows = $db->select_href( ":hashes", "some_field :ORDER BY some_field" );
 @rows = $db->select_href( ":hashes", "f2 f6 f5 :ORDER BY f4" );

 # $rows[0] = [ "key1", { f2 => "val", f6 => "val", f5 => "val" } ]
 # $rows[1] = [ "key2", { f2 => "val", f6 => "val", f5 => "val" } ]
 # $rows[2] = [ "key3", { f2 => "val", f6 => "val", f5 => "val" } ]
 # ...
 # $rows[N] = [ "keyN", { f2 => "val", f6 => "val", f5 => "val" } ]

=head2 select_aref ( ":lists", "select string" )

=head2 select_href ( ":lists", "select string" )

Returns a list containing C<[ key, aref ]> pairs or C<[ key, href ]> pairs
from the hash of lists (HoA).

The C<select_aref> and C<select_href> methods take a select string supporting
field indices and optionally sort modifiers. The syntax for the query string,
between C<:WHERE> and C<:ORDER BY>, is the same as described above.

The modifiers C<:WHERE>, C<:AND>, C<:OR>, C<ORDER BY>, C<ASC>, C<DESC>, and
C<ALPHA> may be mixed case. e.g. C<:Where>

HoA select string:

 "17 15 11 :WHERE 12 > 20 :AND key =~ /foo/ :ORDER BY 10 DESC ALPHA"
 "17 15 11 :WHERE 12 > 40 :AND key =~ /bar/ :ORDER BY key ALPHA"
 "17 15 11 :WHERE 12 > 40 :AND key =~ /bar/"
 "17 15 11"

 * key matches on keys stored in the primary level hash (H)oA
 * above, list indices are given as 17, 15, 11, 12, and 10
 * the shorter form is allowed e.g. "4 > 20 :AND key =~ /baz/"

HoA select synopsis:

 @rows = $db->select_aref( ":lists", "some_index :ORDER BY some_index" );
 @rows = $db->select_aref( ":lists", "2 6 5 :ORDER BY 4" );

 # $rows[0] = [ "key1", [ "2_val", "6_val", "5_val" ] ]
 # $rows[1] = [ "key2", [ "2_val", "6_val", "5_val" ] ]
 # $rows[2] = [ "key3", [ "2_val", "6_val", "5_val" ] ]
 # ...
 # $rows[N] = [ "keyN", [ "2_val", "6_val", "5_val" ] ]

 @rows = $db->select_href( ":lists", "some_index :ORDER BY some_index" );
 @rows = $db->select_href( ":lists", "2 6 5 :ORDER BY 4" );

 # $rows[0] = [ "key1", { 2 => "val", 6 => "val", 5 => "val" } ]
 # $rows[1] = [ "key2", { 2 => "val", 6 => "val", 5 => "val" } ]
 # $rows[2] = [ "key3", { 2 => "val", 6 => "val", 5 => "val" } ]
 # ...
 # $rows[N] = [ "keyN", { 2 => "val", 6 => "val", 5 => "val" } ]

=head1 API DOCUMENTATION - HASHES ( HoH )

=head2 hassign ( key, field, value [, field, value, ... ] )

Clears the hash stored at key, then sets the value of a hash field and returns
its new value. Multiple field_value pairs may be set at once. In that case, the
number of fields is returned. This is equivalent to C<hclear>, C<hset>.

 $val = $db->hassign( "some_key", "field", "value" );
 $len = $db->hassign( "some_key", "f1" => "val1", "f2" => "val2" );

API available since 1.007.

=head2 hclear ( key )

Removes all key-value pairs from the first level hash (H)oH when no arguments
are given. Otherwise, removes all field-value pairs stored at key.

 $db->hclear;
 $db->hclear( "some_key" );

=head2 hdel ( key, field [, field, ... ] )

Deletes one or more hash fields. It returns the value associated with the field
if a single field is given. Otherwise, it returns the number of fields actually
removed from the hash stored at key. A field which does not exist in the hash
is not counted.

 $val = $db->hdel( "some_key", "some_field" );
 $cnt = $db->hdel( "some_key", "field1", "field2" );

=head2 hdel ( key )

Deletes and returns the C<MCE::Shared::Hash> object stored at key or C<undef>
if the key does not exists in the first level hash (H)oH.

 $ha_obj = $db->hdel( "some_key" );

=head2 hexists ( key, field [, field, ... ] )

Determines if a hash field exists. For multiple fields, a truth value is
returned only if all given fields exist in the hash stored at key.

 if ( $db->hexists( "some_key", "some_field" ) ) { ... }
 if ( $db->hexists( "some_key", "f1", "f5" ) ) { ... }

=head2 hexists ( key )

Determines if a key exists in the first level hash (H)oH.

 if ( $db->hexists( "some_key" ) ) { ... }

=head2 hget ( key, field [, field, ... ] )

Gets the values of all given hash fields. The C<undef> value is returned for
fields which do not exists in the hash stored at key. Likewise, the C<undef>
value is returned if the key does not exists in the first level hash (H)oH.

 $val = $db->hget( "some_key", "field" );

 ( $val1, $val2 ) = $db->hget( "some_key", "field1", "field2" );

=head2 hget ( key )

Gets the C<MCE::Shared::Hash> object for the hash stored at key or C<undef> if
the key does not exists in the first level hash (H)oH.

 $ha_obj = $db->hget( "some_key" );

=head2 hkeys ( key, [ field [, field, ... ] ] )

Returns keys stored in the first level hash (H)oH when no arguments are given.
Otherwise, returns the given fields in the hash stored at key. Fields that do
not exist will have the C<undef> value. In scalar context, returns the size of
the object, either the hash stored at key or the first level hash.

 @keys   = $db->hkeys;
 @fields = $db->hkeys( "some_key" );
 @fields = $db->hkeys( "some_key", "field1", "field2" );
 $len    = $db->hkeys( "some_key" );
 $len    = $db->hkeys;

=head2 hkeys ( key, "query string" )

Returns only fields stored at key that match the given criteria. It returns an
empty list if the search found nothing. The syntax for the C<query string> is
described above. In scalar context, returns the size of the resulting list.

 @keys = $db->hkeys( "some_key", "val eq some_value" );
 @keys = $db->hkeys( "some_key", "key eq some_key :AND val =~ /sun|moon/" );
 @keys = $db->hkeys( "some_key", "val eq sun :OR val eq moon );
 $len  = $db->hkeys( "some_key", "key =~ /$pattern/" );

=head2 hkeys ( "query string" )

For the one argument form, the search is applied to keys stored in the primary
hash only (H)oH. Therefore, the C<key> modifier is the only thing possible if
wanting any result.

 @keys = $db->hkeys( "key =~ /$pattern/" );
 $len  = $db->hkeys( "key =~ /$pattern/" );

=head2 hlen ( key [, field ] )

Returns the size of the first level hash (H)oH when no arguments are given.
For the given key, returns the size of hash stored at key or optionally, the
length of the value stored at key-field. It returns the C<undef> value if
either the given key or given field does not exists.

 $len = $db->hlen;
 $len = $db->hlen( $key );
 $len = $db->hlen( $key, $field );

=head2 hpairs ( key, [ field [, field, ... ] ] )

Returns key-value pairs stored in the first level hash (H)oH when no arguments
are given. Otherwise, returns field-value pairs for the given fields in the hash
stored at key. Fields that do not exist will have the C<undef> value. In scalar
context, returns the size of the object, either the hash stored at key or the
first level hash.

 @pairs = $db->hpairs;                 # ( key => href, ... )
 @pairs = $db->hpairs( "some_key" );   # ( field => value, ... )
 @pairs = $db->hpairs( "some_key", "field1", "field2" );
 $len   = $db->hpairs( "some_key" );
 $len   = $db->hpairs;

=head2 hpairs ( key, "query string" )

Returns only field-value pairs stored at key that match the given criteria.
It returns an empty list if the search found nothing. The syntax for the
C<query string> is described above. In scalar context, returns the size of
the resulting list.

 @pairs = $db->hpairs( "some_key", "val eq some_value" );
 @pairs = $db->hpairs( "some_key", "key eq some_key :AND val =~ /sun|moon/" );
 @pairs = $db->hpairs( "some_key", "val eq sun :OR val eq moon" );
 $len   = $db->hpairs( "some_key", "key =~ /$pattern/" );

=head2 hpairs ( "query string" )

For the one argument form, the search is applied to keys stored in the primary
hash only (H)oH. Therefore, the C<key> modifier is the only thing possible if
wanting any result.

 @keys = $db->hpairs( "key =~ /$pattern/" );
 $len  = $db->hpairs( "key =~ /$pattern/" );

=head2 hset ( key, field, value [, field, value, ... ] )

Sets the value of a hash field and returns its new value. Multiple field_value
pairs may be set at once. In that case, the number of fields stored at key is
returned.

 $val = $db->hset( "some_key", "field", "value" );
 $len = $db->hset( "some_key", "f1" => "val1", "f2" => "val2" );

=head2 hsetnx ( key, field, value )

Sets the value of a hash field, only if the field does not exist. Returns a 1
for new field or 0 if the field already exists and no operation was performed.

 $ret = $db->hsetnx( "some_key", "field", "value" );

Current API available since 1.872.

=head2 hshift

Removes and returns the first key-value pair or value in scalar context from
the first-level hash (H)oH. If the C<HASH> is empty, returns the undefined
value.

 ( $key, $href ) = $db->hshift;

 $href = $db->hshift;

=head2 hsort ( "BY key [ ASC | DESC ] [ ALPHA ]" )

=head2 hsort ( "BY field [ ASC | DESC ] [ ALPHA ]" )

Returns sorted keys from the first level hash (H)oH, leaving the elements
intact. In void context, sorts the first level hash in-place. Sorting is
numeric by default.

 @keys = $db->hsort( "BY key" );
 @keys = $db->hsort( "BY some_field" );

 $db->hsort( "BY key" );
 $db->hsort( "BY some_field" );

If the keys or field values contain string values and you want to sort them
lexicographically, specify the C<ALPHA> modifier.

 @keys = $db->hsort( "BY key ALPHA" );
 @keys = $db->hsort( "BY some_field ALPHA" );

 $db->hsort( "BY key ALPHA" );
 $db->hsort( "BY some_field ALPHA" );

The default is C<ASC> for sorting the elements from small to large. In order
to sort from large to small, specify the C<DESC> modifier.

 @keys = $db->hsort( "BY key DESC ALPHA" );
 @keys = $db->hsort( "BY some_field DESC ALPHA" );

 $db->hsort( "BY key DESC ALPHA" );
 $db->hsort( "BY some_field DESC ALPHA" );

=head2 hvals ( key, [ field [, field, ... ] ] )

Returns values stored in the first level hash (H)oH when no arguments are given.
Otherwise, returns values for the given fields in the hash stored at key. Fields
that do not exist will have the C<undef> value. In scalar context, returns the
size of the object, either the hash stored at key or the first level hash.

 @hrefs = $db->hvals;
 @vals  = $db->hvals( "some_key" );
 @vals  = $db->hvals( "some_key", "field1", "field2" );
 $len   = $db->hvals( "some_key" );
 $len   = $db->hvals;

=head2 hvals ( key, "query string" )

Returns only values stored at key that match the given criteria. It returns an
empty list if the search found nothing. The syntax for the C<query string> is
described above. In scalar context, returns the size of the resulting list.

 @vals = $db->hvals( "some_key", "val eq some_value" );
 @vals = $db->hvals( "some_key", "key eq some_key :AND val =~ /sun|moon/" );
 @vals = $db->hvals( "some_key", "val eq sun :OR val eq moon" );
 $len  = $db->hvals( "some_key", "key =~ /$pattern/" );

=head2 hvals ( "query string" )

For the one argument form, the search is applied to keys stored in the primary
hash only (H)oH. Therefore, the C<key> modifier is the only thing possible if
wanting any result.

 @keys = $db->hvals( "key =~ /$pattern/" );
 $len  = $db->hvals( "key =~ /$pattern/" );

=head1 SUGAR METHODS - HASHES ( HoH )

This module is equipped with sugar methods to not have to call C<set>
and C<get> explicitly. In shared context, the benefit is atomicity and
reduction in inter-process communication.

=head2 happend ( key, field, string )

Appends a value to key-field and returns its new length.

 $len = $db->happend( $key, $field, "foo" );

=head2 hdecr ( key, field )

Decrements the value of key-field by one and returns its new value.

 $num = $db->hdecr( $key, $field );

=head2 hdecrby ( key, field, number )

Decrements the value of key-field by the given number and returns its new value.

 $num = $db->hdecrby( $key, $field, 2 );

=head2 hgetdecr ( key, field )

Decrements the value of key-field by one and returns its old value.

 $old = $db->hgetdecr( $key, $field );

=head2 hgetincr ( key, field )

Increments the value of key-field by one and returns its old value.

 $old = $db->hgetincr( $key, $field );

=head2 hgetset ( key, field, value )

Sets the value of key-field and returns its old value.

 $old = $db->hgetset( $key, $field, "baz" );

=head2 hincr ( key, field )

Increments the value of key-field by one and returns its new value.

 $num = $db->hincr( $key, $field );

=head2 hincrby ( key, field, number )

Increments the value of key-field by the given number and returns its new value.

 $num = $db->hincrby( $key, $field, 2 );

=head1 API DOCUMENTATION - LISTS ( HoA )

=head2 lassign ( key, value [, value, ... ] )

Clears the list stored at key, then prepends one or multiple values and returns
the new length. This is equivalent to C<lclear>, C<lpush>.

 $len = $db->lassign( "some_key", "val1", "val2" );

API available since 1.007.

=head2 lclear ( key )

Removes all key-value pairs from the first level hash (H)oA when no arguments
are given. Otherwise, removes all elements from the list stored at key.

 $db->lclear;
 $db->lclear( "some_key" );

=head2 ldel ( key, index [, index, ... ] )

Deletes one or more elements by their indices. It returns the value associated
with the index if a single index is given. Otherwise, it returns the number of
elements actually removed from the list stored at key. An index which does not
exists in the list is not counted.

 $val = $db->ldel( "some_key", 20 );
 $cnt = $db->ldel( "some_key", 0, 1 );

=head2 ldel ( key )

Deletes and returns the C<MCE::Shared::Array> object stored at key or C<undef>
if the key does not exists in the first level hash (H)oA.

 $ar_obj = $db->ldel( "some_key" );

=head2 lexists ( key, index [, index, ... ] )

Determines if elements by their indices exist in the list. For multiple indices,
a truth value is returned only if all given indices exist in the list stored at
key. The behavior is strongly tied to the use of delete on lists.

 $db->lset( "some_key", 0, "value0" );
 $db->lset( "some_key", 1, "value1" );
 $db->lset( "some_key", 2, "value2" );
 $db->lset( "some_key", 3, "value3" );

 $db->lexists( "some_key", 2 );     # True
 $db->lexists( "some_key", 2, 3 );  # True
 $db->ldel   ( "some_key", 2 );     # value2

 $db->lexists( "some_key", 2 );     # False
 $db->lexists( "some_key", 2, 3 );  # False
 $db->lexists( "some_key", 3 );     # True

=head2 lexists ( key )

Determines if a key exists in the first level hash (H)oA.

 if ( $db->lexists( "some_key" ) ) { ... }

=head2 lget ( key, index [, index, ... ] )

Gets the values of all given list indices. The C<undef> value is returned for
indices which do not exists in the list stored at key. Likewise, the C<undef>
value is returned if the key does not exists in the first level hash (H)oA.

 $val = $db->lget( "some_key", 20 );

 ( $val1, $val2 ) = $db->lget( "some_key", 0, 1 );

=head2 lget ( key )

Gets the C<MCE::Shared::Array> object for the list stored at key or C<undef> if
the key does not exists in the first level hash (H)oA.

 $ar_obj = $db->lget( "some_key" );

=head2 lkeys ( key, [ index [, index, ... ] ] )

Returns keys stored in the first level hash (H)oA when no arguments are given.
Otherwise, returns the given indices in the list stored at key. Indices that do
not exist will have the C<undef> value. In scalar context, returns the size of
the object, either the list stored at key or the first level hash.

 @keys    = $db->lkeys;
 @indices = $db->lkeys( "some_key" );
 @indices = $db->lkeys( "some_key", 0, 1 );
 $len     = $db->lkeys( "some_key" );
 $len     = $db->lkeys;

=head2 lkeys ( key, "query string" )

Returns only indices stored at key that match the given criteria. It returns an
empty list if the search found nothing. The syntax for the C<query string> is
described above. In scalar context, returns the size of the resulting list.

 @keys = $db->lkeys( "some_key", "val eq some_value" );
 @keys = $db->lkeys( "some_key", "key >= 50 :AND val =~ /sun|moon/" );
 @keys = $db->lkeys( "some_key", "val eq sun :OR val eq moon" );
 $len  = $db->lkeys( "some_key", "key =~ /$pattern/" );

=head2 lkeys ( "query string" )

For the one argument form, the search is applied to keys stored in the primary
hash only (H)oA. Therefore, the C<key> modifier is the only thing possible if
wanting any result.

 @keys = $db->lkeys( "key =~ /$pattern/" );
 $len  = $db->lkeys( "key =~ /$pattern/" );

=head2 llen ( key [, index ] )

Returns the size of the first level hash (H)oA when no arguments are given.
For the given index, returns the size of list stored at key or optionally,
the length of the value stored at key-index. It returns the C<undef> value
if either the given key or given index does not exists.

 $len = $db->llen;
 $len = $db->llen( $key );
 $len = $db->llen( $key, 0 );

=head2 lpairs ( key, [ index [, index, ... ] ] )

Returns key-value pairs stored in the first level hash (H)oA when no arguments
are given. Otherwise, returns index-value pairs for the given indices in the
list stored at key. Indices that do not exist will have the C<undef> value.
In scalar context, returns the size of the object, either the list stored
at key or the first level hash.

 @pairs = $db->lpairs;                 # ( key => aref, ... )
 @pairs = $db->lpairs( "some_key" );   # ( index => value, ... )
 @pairs = $db->lpairs( "some_key", 0, 1 );
 $len   = $db->lpairs( "some_key" );
 $len   = $db->lpairs;

=head2 lpairs ( key, "query string" )

Returns only index-value pairs stored at key that match the given criteria.
It returns an empty list if the search found nothing. The syntax for the
C<query string> is described above. In scalar context, returns the size of
the resulting list.

 @pairs = $db->lpairs( "some_key", "val eq some_value" );
 @pairs = $db->lpairs( "some_key", "key >= 50 :AND val =~ /sun|moon/" );
 @pairs = $db->lpairs( "some_key", "val eq sun :OR val eq moon" );
 $len   = $db->lpairs( "some_key", "key =~ /$pattern/" );

=head2 lpairs ( "query string" )

For the one argument form, the search is applied to keys stored in the primary
hash only (H)oA. Therefore, the C<key> modifier is the only thing possible if
wanting any result.

 @keys = $db->lpairs( "key =~ /$pattern/" );
 $len  = $db->lpairs( "key =~ /$pattern/" );

=head2 lpop ( key )

Removes and returns the first value of the list stored at key. If there are
no elements in the list, returns the undefined value.

 $val = $db->lpop( $key );

=head2 lpush ( key, value [, value, ... ] )

Prepends one or multiple values to the head of the list stored at key and
returns the new length.

 $len = $db->lpush( "some_key", "val1", "val2" );

=head2 lrange ( key, start, stop )

Returns the specified elements of the list stored at key. The offsets C<start>
and C<stop> can also be negative numbers indicating offsets starting at the
end of the list.

An empty list is returned if C<start> is larger than the end of the list.
C<stop> is set to the last index of the list if larger than the actual end
of the list.

 @list = $db->lrange( "some_key", 20, 29 );
 @list = $db->lrange( "some_key", -4, -1 );

=head2 lset ( key, index, value [, index, value, ... ] )

Sets the value of an element in a list by its index and returns its new value.
Multiple index_value pairs may be set all at once. In that case, the length of
the list is returned.

 $val = $db->lset( "some_key", 2, "value" );
 $len = $db->lset( "some_key", 0 => "val1", 1 => "val2" );

=head2 lshift

Removes and returns the first key-value pair or value in scalar context from
the first-level hash (H)oA. If the C<HASH> is empty, returns the undefined
value. See C<lpop> to shift the first value of the list stored at key.

 ( $key, $aref ) = $db->lshift;

 $aref = $db->lshift;

=head2 lsort ( "BY key [ ASC | DESC ] [ ALPHA ]" )

=head2 lsort ( "BY index [ ASC | DESC ] [ ALPHA ]" )

Returns sorted keys from the first level hash (H)oA, leaving the elements
intact. In void context, sorts the first level hash in-place. Sorting is
numeric by default.

 @keys = $db->lsort( "BY key" );
 @keys = $db->lsort( "BY some_index" );

 $db->lsort( "BY key" );
 $db->lsort( "BY some_index" );

If the keys or field values given by index contain string values and you want
to sort them lexicographically, specify the C<ALPHA> modifier.

 @keys = $db->lsort( "BY key ALPHA" );
 @keys = $db->lsort( "BY some_index ALPHA" );

 $db->lsort( "BY key ALPHA" );
 $db->lsort( "BY some_index ALPHA" );

The default is C<ASC> for sorting the elements from small to large. In order
to sort from large to small, specify the C<DESC> modifier.

 @keys = $db->lsort( "BY key DESC ALPHA" );
 @keys = $db->lsort( "BY some_index DESC ALPHA" );

 $db->lsort( "BY key DESC ALPHA" );
 $db->lsort( "BY some_index DESC ALPHA" );

=head2 lsort ( key, "BY key [ ASC | DESC ] [ ALPHA ]" )

=head2 lsort ( key, "BY val [ ASC | DESC ] [ ALPHA ]" )

The two argument form has similar functionality. Here, sorting is applied to
the list stored at key. For example, the following reverses the list stored
at the given key. The C<BY key> modifier refers to the indices in the list
and not the values.

 $db->lsort( "some_key", "BY key DESC" );

=head2 lsplice ( key, offset [, length [, list ] ] )

Removes the elements designated by C<offset> and C<length> from the array
stored at key, and replaces them with the elements of C<list>, if any.
The behavior is similar to the Perl C<splice> function.

 @items = $db->lsplice( "some_key", 20, 2, @list );
 @items = $db->lsplice( "some_key", 20, 2 );
 @items = $db->lsplice( "some_key", 20 );

=head2 lvals ( key, [ index [, index, ... ] ] )

Returns values stored in the first level hash (H)oA when no arguments are
given. Otherwise, returns values for the given indices in the list stored
at key. Indices that do not exist will have the C<undef> value. In scalar
context, returns the size of the object, either the list stored at key or
the first level hash.

 @arefs = $db->lvals;
 @vals  = $db->lvals( "some_key" );
 @vals  = $db->lvals( "some_key", 0, 1 );
 $len   = $db->lvals( "some_key" );
 $len   = $db->lvals;

=head2 lvals ( key, "query string" )

Returns only values stored at key that match the given criteria. It returns an
empty list if the search found nothing. The syntax for the C<query string> is
described above. In scalar context, returns the size of the resulting list.

 @keys = $db->lvals( "some_key", "val eq some_value" );
 @keys = $db->lvals( "some_key", "key >= 50 :AND val =~ /sun|moon/" );
 @keys = $db->lvals( "some_key", "val eq sun :OR val eq moon" );
 $len  = $db->lvals( "some_key", "key =~ /$pattern/" );

=head2 lvals ( "query string" )

For the one argument form, the search is applied to keys stored in the primary
hash only (H)oA. Therefore, the C<key> modifier is the only thing possible if
wanting any result.

 @keys = $db->lvals( "key =~ /$pattern/" );
 $len  = $db->lvals( "key =~ /$pattern/" );

=head2 rpop ( key )

Removes and returns the last value of the list stored at key. If there are
no elements in the list, returns the undefined value.

 $val = $db->rpop( $key );

=head2 rpush ( key, value [, value, ... ] )

Appends one or multiple values to the tail of the list stored at key and
returns the new length.

 $len = $db->rpush( "some_key", "val1", "val2" );

=head1 SUGAR METHODS - LISTS ( HoA )

This module is equipped with sugar methods to not have to call C<set>
and C<get> explicitly. In shared context, the benefit is atomicity and
reduction in inter-process communication.

=head2 lappend ( key, index, string )

Appends a value to key-index and returns its new length.

 $len = $db->lappend( $key, 0, "foo" );

=head2 ldecr ( key, index )

Decrements the value of key-index by one and returns its new value.

 $num = $db->ldecr( $key, 0 );

=head2 ldecrby ( key, index, number )

Decrements the value of key-index by the given number and returns its new value.

 $num = $db->ldecrby( $key, 0, 2 );

=head2 lgetdecr ( key, index )

Decrements the value of key-index by one and returns its old value.

 $old = $db->lgetdecr( $key, 0 );

=head2 lgetincr ( key, index )

Increments the value of key-index by one and returns its old value.

 $old = $db->lgetincr( $key, 0 );

=head2 lgetset ( key, index, value )

Sets the value of key-index and returns its old value.

 $old = $db->lgetset( $key, 0, "baz" );

=head2 lincr ( key, index )

Increments the value of key-index by one and returns its new value.

 $num = $db->lincr( $key, 0 );

=head2 lincrby ( key, index, number )

Increments the value of key-index by the given number and returns its new value.

 $num = $db->lincrby( $key, 0, 2 );

=head1 CREDITS

The implementation is inspired by various Redis Hash/List primitives at
L<https://redis.io/commands>.

=head1 INDEX

L<MCE|MCE>, L<MCE::Hobo>, L<MCE::Shared>

=head1 AUTHOR

Mario E. Roy, S<E<lt>marioeroy AT gmail DOT comE<gt>>

=cut

