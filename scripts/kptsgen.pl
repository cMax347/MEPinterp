#!/usr/bin/perl -w
#
# Generates a kpts file for band-structure calculations with Fleur. The lattice
# names are the same as in the input generator, and the high-symmetry k-points
# follow standard notation. The command-line arguments -l and -k are mandatory.
#
# Author: A. Schindlmayr (A.Schindlmayr@fz-juelich.de)
# In order to include a new lattice:
# - Add the three-letter short name of the lattice to @lattices.
# - Define high-symmetry points in subroutine highsymmetry.

use strict;
use Getopt::Std;
use vars qw ($opt_l $opt_k $opt_o $opt_h);

# List of lattices and synonyms
my @lattices = ( "cub", "fcc", "bcc", "hcp", "tet", "orP" );
my %synonyms = ( cub => "simple-cubic sc cP cubic-P",
                 fcc => "face-centered-cubic cF cubic-F",
                 bcc => "body-centered-cubic cI cubic-I",
                 hcp => "hexagonal hP hexagonal-P",
                 rho => "rhombohedral hr r R hexagonal-R trigonal",
                 tet => "simple-tetragonal st tP tetragonal-P",
                 bct => "body-centered-tetragonal tI tetragonal-I",
                 orP => "simple-orthorhombic oP orthorhombic-P",
                 orF => "face-centered-orthorhombic oF orthorhombic-F",
                 orI => "body-centered-orthorhombic oI orthorhombic-I",
                 orC => "base-centered-orthorhombic oC oS orthorhombic-C orthorhombic-S",
                 moP => "simple-monoclinic mP monoclinic-P",
                 moC => "centered-monoclinic mC monoclinic-C",
                 tcl => "triclinic aP",
                 hdp => "hexagonal2 hexagonal-2",
                 trg => "rhombohedral2 hR2 r2 R2 hexagonal-R2",
                 orA => "base-centered-orthorhombic2 oA orthorhombic-A",
                 orB => "base-centered-orthorhombic3 oB orthorhombic-B",
                 moA => "centered-monoclinic2 mA monoclinic-A",
                 moB => "centered-monoclinic3 mB monoclinic-B"
               );

# Get options from command line
$opt_o = "kpts";
getopts('l:k:o:h');

# Print help message if -h
if ($opt_h) {
  print <<USAGE_END;

Generates a kpts file for band-structure calculations with Fleur. The lattice
names are the same as in the input generator, and the high-symmetry k-points
follow standard notation. The command-line arguments -l and -k are mandatory.

Usage: $0 -l <lattice> -k <path> [-o <outfile>] [-h]

Options:
  -l   lattice, e.g. face-centered-cubic or fcc
  -k   string with k-points and number of steps, e.g. "L 13 Gamma 15 X"
  -o   output file (default: kpts)
  -h   help

The following lattices and high-symmetry points are currently implemented.
USAGE_END
  foreach $opt_l (@lattices) {
    my %k = &highsymmetry;
    print "  ",$opt_l,(exists $synonyms{$opt_l})?" ($synonyms{$opt_l})":"",
      ": ",join (", ",sort (keys %k)),"\n";
  }
  exit 0;
}

# Check whether a valid lattice is specified with -l
die "No lattice specified. Use -h for help.\n" unless ($opt_l);
my %k = &highsymmetry;

# Check whether a valid path is specified with -k
die "No path specified. Use -h for help.\n" unless ($opt_k);
my ($kpoint1,@path) = split (" ",$opt_k);
die "High-symmetry k-point $kpoint1 not defined in $opt_l lattice.\n"
  unless (exists $k{$kpoint1});
my @output = ( sprintf ("%10.5f%10.5f%10.5f%10.5f   $kpoint1\n",
  ${$k{$kpoint1}}[0],${$k{$kpoint1}}[1],${$k{$kpoint1}}[2],1.0) );
while (@path) {
  my $nsteps = shift (@path);
  die "Number of steps must be positive.\n" unless ($nsteps>0);
  die "No terminating k-point specified.\n" unless (@path);
  my $kpoint0 = $kpoint1; $kpoint1 = shift (@path);
  die "High-symmetry k-point $kpoint1 not defined in $opt_l lattice.\n"
    unless (exists $k{$kpoint1});
  for (my $istep = 1; $istep <= $nsteps; $istep++) {
    push (@output, sprintf ("%10.5f%10.5f%10.5f%10.5f%s\n",
      ${$k{$kpoint0}}[0]+($istep/$nsteps)*(${$k{$kpoint1}}[0]-${$k{$kpoint0}}[0]),
      ${$k{$kpoint0}}[1]+($istep/$nsteps)*(${$k{$kpoint1}}[1]-${$k{$kpoint0}}[1]),
      ${$k{$kpoint0}}[2]+($istep/$nsteps)*(${$k{$kpoint1}}[2]-${$k{$kpoint0}}[2]),
      1.0, ($istep==$nsteps)?"   $kpoint1":"" ));
  }
}

# Write kpts file
open (OUT,">$opt_o") or die "Output file $opt_o cannot be opened.\n";
printf OUT ("%5d%20.10f\n",$#output+1,1.0);
print OUT @output;
close (OUT);


# Define high-symmetry points
sub highsymmetry {
  my $lattice = $opt_l;
  foreach $_ (keys %synonyms) {
    @_ = split (" ",$synonyms{$_});
    while (@_) {
      $lattice = $_ if ($opt_l eq shift (@_));
    }
  }
  if ($lattice eq "cub") {            # simple-cubic
    %k = (
           Gamma => [  0 ,  0 ,  0 ],
           M     => [ 1/2, 1/2,  0 ],
           R     => [ 1/2, 1/2, 1/2],
           X     => [ 1/2,  0 ,  0 ]
         );
  } elsif ($lattice eq "fcc") {       # face-centered-cubic
    %k = (
           Gamma => [  0 ,  0 ,  0 ],
           K     => [ 3/8, 3/8, 3/4],
           L     => [ 1/2, 1/2, 1/2],
           U     => [ 1/4, 5/8, 5/8],
           W     => [ 1/4, 1/2, 3/4],
           X     => [  0 , 1/2, 1/2],
         );
  } elsif ($lattice eq "bcc") {       # body-centered-cubic
    %k = (
           Gamma => [  0 ,  0 ,  0 ],
           H     => [-1/2, 1/2, 1/2],
           N     => [  0 ,  0 , 1/2],
           P     => [ 1/4, 1/4, 1/4]
         );
  } elsif ($lattice eq "hcp") {       # hexagonal
    %k = (
           A     => [  0 ,  0 , 1/2],
           Gamma => [  0 ,  0 ,  0 ],
           H     => [ 1/3, 1/3, 1/2],
           K     => [ 1/3, 1/3,  0 ],
           L     => [  0 , 1/2, 1/2],
           M     => [  0 , 1/2,  0 ]
         );
  } elsif ($lattice eq "tet") {       # simple-tetragonal
    %k = (
           A     => [ 1/2, 1/2, 1/2],
           Gamma => [  0 ,  0 ,  0 ],
           M     => [ 1/2, 1/2,  0 ],
           R     => [ 1/2,  0 , 1/2],
           X     => [ 1/2,  0 ,  0 ],
           Z     => [  0 ,  0 , 1/2]
         );
  } elsif ($lattice eq "orP") {       # simple-orthorhombic
    %k = (
           Gamma => [  0 ,  0 ,  0 ],
           R     => [ 1/2, 1/2, 1/2],
           S     => [ 1/2, 1/2,  0 ],
           T     => [  0 , 1/2, 1/2],
           U     => [ 1/2,  0 , 1/2],
           X     => [ 1/2,  0 ,  0 ],
           Y     => [  0 , 1/2,  0 ],
           Z     => [  0 ,  0 , 1/2]
         );
    } elsif ($lattice eq "jpH") {       # JANS 3Q MODEL
    %k = (
           Gamma => [  0 ,  0 ,  0 ],
           R     => [  0 ,  0 , 1/2],
           X     => [ 1/2,  0 ,  0 ],
           M     => [ 1/2, 1/3,  0 ]
         );
  } else {
    die "Unknown lattice: $opt_l.\n";
  }
}
