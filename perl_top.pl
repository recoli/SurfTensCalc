=comment
 *  This file is part of the SurfTensCalc program.
 *  Copyright (C) 2012 Xin Li

 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.

 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.

 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
=cut

#!/usr/bin/perl 

   use strict;
   use warnings;

#  printf title
   printf "\n";
   printf "***************************************************************\n";
   printf "*      Perl script to convert gmx top file to param.txt       *\n";
   printf "*                        by  Xin Li                           *\n";
   printf "*         TheoChem & Biology, KTH, Stockholm, Sweden          *\n";
   printf "***************************************************************\n";
   printf "\n";
   printf "                     Version 2012-03-16                        \n";
   printf "\n";
   printf "Usage: perl perl_top.pl -top [gmx_top_file] -trr [gmx_trr_file]\n";
   printf "\n";

#  options
#  Note: change array to hash
   my %options = @ARGV;

   die "Error: undefined -top option!\n" unless (defined $options{"-top"});
   my $gmx_top = $options{"-top"};
   die "Error: cannot find file $gmx_top !\n" unless (-e $gmx_top);

   die "Error: undefined -trr option!\n" unless (defined $options{"-trr"});
   my $gmx_trr = $options{"-trr"};
   die "Error: cannot find file $gmx_trr !\n" unless (-e $gmx_trr);

   printf "GMX_top_file   $gmx_top\n";
   printf "GMX_trr_file   $gmx_trr\n";
   printf "\n";

#  check:
#  OPLS, AMBER03 force field
#  fudgeLJ, fudgeQQ
#  SPC/E water model
   my $use_opls = 0; my $use_amber03 = 0;
   my $fudge_LJ; my $fudge_QQ; my $comb_rule;
   my $use_spce = 0; 

   open TOP,"$gmx_top" 
      or die "Error: cannot open file $gmx_top !\n";
   while ( <TOP> ) 
   {
      if ( /^\#include\s+\"oplsaa\.ff\/forcefield\.itp\"/ 
           or /^\#include\s+\"ffoplsaa\.itp\"/ )
      {
         $use_opls = 1;
         $fudge_LJ = 0.5;
         $fudge_QQ = 0.5;
         $comb_rule = 3;
      }
      if ( /^\#include\s+\"amber03\.ff\/forcefield\.itp\"/ 
           or /^\#include\s+\"ffamber03\.itp\"/ )
      {
         $use_amber03 = 1;
         $fudge_LJ = 0.5;
         $fudge_QQ = 1.0/1.2;
         $comb_rule = 2;
      }

      if ( /^\#include\s+\"oplsaa\.ff\/spce\.itp\"/ 
           or /^\#include\s+\"amber03\.ff\/spce\.itp\"/
           or /^\#include\s+\"spce\.itp\"/)
      {
         $use_spce = 1;
      }

      if ( /^\[\s*defaults\s*\]/ )
      {
         while ( <TOP> )
         {
#           skip comments, # lines and blank lines
            next if ( /^;/ or /^\#/ or /^\s*$/ );
            next if ( /^\[\s*defaults\s*\]/ );
#           do not read next block
            last if( /^\[/ );
#           read fudgeLJ and fudgeQQ
            @_ = (split);
            $comb_rule = $_[1];
            if ( defined $_[3] and $_[4] )
            {
               $fudge_LJ = $_[3];
               $fudge_QQ = $_[4];
            }
            else
            {
               $fudge_LJ = 1.0;
               $fudge_QQ = 1.0;
            }
         }
      }

   }
   close TOP;

#  save atomtypes in another file
   open TOP,"$gmx_top" 
      or die "Error: cannot open file $gmx_top !\n";
   open ITP,">$gmx_top.atomtypes.itp"
      or die "Error: cannot creat file $gmx_top.atomtypes.itp !\n";
   while ( <TOP> ) 
   {
      if ( /^\[\s*atomtypes\s*\]/ )
      {
         printf ITP;
         while ( <TOP> )
         {
#           skip comments, # lines and blank lines
            next if ( /^;/ or /^\#/ or /^\s*$/ );
            next if ( /^\[\s*atomtypes\s*\]/ );
#           do not read next block
            last if( /^\[/ );
#           write to file
            printf ITP;
         }
      }
   }
   close ITP;
   close TOP;

#  check force field itp file
   if ( $use_opls == 1 )
   {
      die "Error: cannot find file ffoplsaanb.itp!\n"
         unless (-e "itp_files/ffoplsaanb.itp") ;
      die "Error: cannot find file ffoplsaabd.itp!\n"
         unless (-e "itp_files/ffoplsaabd.itp") ;
      printf "Force_Field   OPLS-AA\nffoplsaanb.itp\nffoplsaabd.itp\n";
   }
   if ( $use_amber03 == 1 )
   {
      die "Error: cannot find file ffamber03nb.itp!\n"
         unless (-e "itp_files/ffamber03nb.itp") ;
      die "Error: cannot find file ffamber03bd.itp!\n"
         unless (-e "itp_files/ffamber03bd.itp") ;
      printf "Force_Field   Amber03\nffamber03nb.itp\nffamber03bd.itp\n";
   }
   printf "\n";

   if ( $use_spce == 1 )
   {
      printf "Water_Model   SPC/E\n";
   }
   printf "\n";

   printf "Combination_Rule   $comb_rule\n";
   if ( $comb_rule == 2 )
   {
      printf "SIG_ij = 0.5*( SIG_i + SIG_j )\n";
      printf "EPS_ij = sqrt( EPS_i * EPS_j )\n";
   }
   elsif ( $comb_rule == 3 )
   {
      printf "SIG_ij = sqrt( SIG_i * SIG_j )\n";
      printf "EPS_ij = sqrt( EPS_i * EPS_j )\n";
   }
   else
   {
      die "Error: unrecognized comb_rule: $comb_rule !\n";
   }

   printf "\n";

#  open file for writing
   open PAR,">param.txt" or die;

   my $tmp_string;
   my $mol_info;  $mol_info = "";
   my $atom_info;  $atom_info = "";

   my $date = `date`;
   chomp ($date);
   printf PAR "%-15s%-10d\n", "Comb_Rule", $comb_rule;
   printf PAR "%-15s%-15.8f%-15.8f\n", "fudge_LJ_QQ", $fudge_LJ, $fudge_QQ;
   printf PAR "%-15s%-s\n", "trr_file", $gmx_trr;

#  read and write molecular information
   my $mol; my $nMols;
   my @molName; my @molNumber;

   $mol = 0;
   open TOP,"$gmx_top" 
      or die "Error: cannot open gmx top file $gmx_top !\n";
   while ( <TOP> )
   {
      if ( /^\[\s*molecules\s*\]/ )
      {
         while ( <TOP> )
         {
#           skip comments, # lines and blank lines
            next if ( /^;/ or /^\#/ or /^\s*$/ );
            next if ( /^\[\s*molecules\s*\]/ );
#           do not read next block
            last if( /^\[/ );
#           read molecular name and number
            $mol++;
            $molName[$mol] = (split)[0];
            $molName[$mol] = substr $molName[$mol],0,8;
            $molNumber[$mol] = (split)[1];
         }
      }
   }
   $nMols = $mol;

   close TOP;

   $tmp_string = sprintf "%-15s%-10d\n", "MOLECULES", $nMols;
   $mol_info .= $tmp_string;

#  write molecule info on screen
   printf "Molecule_Types   $nMols\n";
   for $mol (1..$nMols)
   {
      printf "%s ", $molName[$mol];
   }
   printf "\n\n";


#  read and write atomic information
   open TOP,"$gmx_top" 
      or die "Error: cannot open gmx top file $gmx_top !\n";

#  variables and arrays
   my ($iType, $jType, $nTypes);
   my ($isNewType, $thisType);

   my ( $atom, $nAtoms );

   my ( @atomResNr, @atomResName );
   my ( @atomName, @atomType, @typeID, @atomCharge, @atomMass );

   my ( $atom_i, $atom_j, $atom_k, $atom_l, $funct );
   my ( $ff_atom_i, $ff_atom_j, $ff_atom_k, $ff_atom_l, $ff_funct );
   my ( $b0, $kb, $D, $beta, $pair, $a0, $ka, $dih );
   my ( @bond_info, @pair_info, @angle_info, @dihedral_info );
   my ( $bond_found, $angle_found, $dihedral_found );
   my ( @nBonds, @nPairs, @nAngles, @nDihedrals );

   my ( $potential, $ii, @coeff );

   my (@typeName, @typeSIG, @typeEPS);

   my $typeFound;

#  initialize number of types and number of molecules
   $mol = 0;
   $nTypes = 0;

   while ( <TOP> )
   {

#     [ atoms ] block: nonbonded parameters
      if ( /^\[\s*atoms\s*\]/ )
      {
         $mol++;

         $nBonds[$mol] = 0;
         $bond_info[$mol] = "";

         $nPairs[$mol] = 0;
         $pair_info[$mol] = "";

         $nAngles[$mol] = 0;
         $angle_info[$mol] = "";

         $nDihedrals[$mol] = 0;
         $dihedral_info[$mol] = "";

         while ( <TOP> )
         {
#           skip comments, # line or blank lines
            next if ( /^;/ or /^\#/ or /^\s*$/ );
            next if ( /^\[\s*atoms\s*\]/ );
#           do not read next block
            last if( /^\[/ );
            $_ = (split ";",$_)[0];
            @_ = (split);

#           read atom ID, resnr, resname, name, charge, mass
            $atom = $_[0];
            $nAtoms = $atom;
            $atomType[$atom]    = $_[1];
            $atomResNr[$atom]   = $_[2];
            $atomResName[$atom] = $_[3];
            $atomName[$atom]    = $_[4];
            $atomCharge[$atom]  = $_[6];
            $atomMass[$atom]    = $_[7];

#           check if this type is a new type 
            $isNewType = 1;
            for $iType (1..$nTypes) 
            {
#              NOTE: use "eq" to compare two strings (not "==")!
               if ( $typeName[$iType] eq $_[1] ) 
               {
                  $isNewType = 0;
                  $thisType = $iType;
                  last;
               }
            }

#           add a new type if this is new
            if ( $isNewType == 1 ) 
            {
               $nTypes++;
               $typeName[$nTypes] = $_[1];

               $typeFound = 0;

#              search this new type in the itp files
               open ITP,"$gmx_top.atomtypes.itp"
                  or die "Error: cannot open itp file $gmx_top.atomtypes.itp !\n";
               while ( <ITP> )
               {
                  $_ = (split ";",$_)[0];
                  @_ = (split);
                  if ( /^\s*$typeName[$nTypes]\s+/ ) 
                  {
                     $typeSIG[$nTypes] = $_[@_-2];
                     $typeEPS[$nTypes] = $_[@_-1];
                     $typeFound = 1;
                  }
               }
               close ITP;

               if ( $use_opls==1 )
               {
                  open ITP,"itp_files/ffoplsaanb.itp"
                     or die "Error: cannot open itp file ffoplsaanb.itp !\n";
                  while ( <ITP> )
                  {
                     $_ = (split ";",$_)[0];
                     @_ = (split);
                     if ( /^\s*$typeName[$nTypes]\s+/ ) 
                     {
                        $typeSIG[$nTypes] = $_[@_-2];
                        $typeEPS[$nTypes] = $_[@_-1];
                        $typeFound = 1;
                     }
                  }
                  close ITP;
               }

               if ( $use_amber03==1 )
               {
                  open ITP,"itp_files/ffamber03nb.itp"
                     or die "Error: cannot open itp file ffamber03nb.itp !\n";
                  while ( <ITP> )
                  {  
                     $_ = (split ";",$_)[0];
                     @_ = (split);
                     if ( /^\s*$typeName[$nTypes]\s+/ ) 
                     {
                        $typeSIG[$nTypes] = $_[@_-2];
                        $typeEPS[$nTypes] = $_[@_-1];
                        $typeFound = 1;
                     }
                  }
                  close ITP;
               }

               if ( $typeFound == 0 )
               {
                  die "Error: cannot find atom type $typeName[$nTypes] in itp file!\n";
               }
               else
               {
#                 get the index of this type
                  $thisType = $nTypes;
               }
            }

#           assign type ID to atom
            $typeID[$atom] = $thisType;

#           change OPLS atomtype to bondtype
            if ( $use_opls==1 )
            {
               open ITP,"itp_files/ffoplsaanb.itp"
                  or die "Error: cannot open itp file ffoplsaanb.itp !\n";
               while ( <ITP> )
               {
                  next if ( /^;/ or /^\[/ or /^\#/ or /^\s*$/ );
                  $_ = (split ";",$_)[0];
                  @_ = (split);
                  $atomType[$atom] = $_[1] if $atomType[$atom] eq $_[0];
               }
               close ITP;
            }

         }

#        write information of a molecule
#        Note: $typeID[$atom]-1 is printed so that it starts from 0.
#              This is beneficial for C program to read.

         $tmp_string = sprintf "%-15s%-10d%-10d\n", 
                       $molName[$mol], $molNumber[$mol], $atomResNr[$nAtoms];
         $mol_info .= $tmp_string;

         $tmp_string = sprintf "%-20s%-10d\n", $molName[$mol]."_ATOMS", $nAtoms;
         $atom_info .= $tmp_string;
         for $atom (1..$nAtoms)
         {
            $tmp_string = sprintf "%13.6f%13.6f%15.8f%15.8f%8d%8d\n",
                          $atomCharge[$atom], $atomMass[$atom],
                          $typeSIG[$typeID[$atom]], $typeEPS[$typeID[$atom]],
                          $typeID[$atom]-1, $atomResNr[$atom]-1;
            $atom_info .= $tmp_string;
         }

#        do not read bonded info if there is only one residue in the molecule
         #next if $atomResNr[$nAtoms]==1;

      }


#     [ bonds ] block: bond stretching
      if ( /^\[\s*bonds\s*\]/ )
      {
         while ( <TOP> )
         {
#           skip comments, # lines and blank lines
            next if ( /^;/ or /^\#/ or /^\s*$/ );
            next if ( /^\[\s*bonds\s*\]/ );
#           do not read next block
            last if( /^\[/ );
#           read bond i, j, funct
            $_ = (split ";",$_)[0];
            @_ = (split);
            $atom_i = $_[0];
            $atom_j = $_[1];
            $funct = $_[2];

            #next if $atomResNr[$atom_i]==$atomResNr[$atom_j];

            if ( $funct==1 or $funct==3 or $funct==5 )
            {
               $nBonds[$mol]++;
            }
            else 
            {
               die "Error: unsupported bond function type: $funct !\n";
            }
#
#           Bond for exclusion
#
            if ( $funct == 5 )
            {
               $tmp_string = sprintf "%d %d   %d %d   %d ;  %s-%s  %s-%s\n",
                                     $atom_i-1, $atom_j-1,  
                                     $atomResNr[$atom_i]-1, $atomResNr[$atom_j]-1,
                                     $funct, 
                                     $atomName[$atom_i], $atomName[$atom_j],
                                     $atomType[$atom_i], $atomType[$atom_j];
               $bond_info[$mol] .= $tmp_string;
            }
#
#           Morse potential
#
            if ( $funct == 3 )
            {
               if ( defined $_[3] and defined $_[4] and defined $_[5] )
               {
                  $b0 = $_[3];
                  $D  = $_[4];
                  $beta = $_[5];
               }
               else
               {
                  printf "Error: bad bond parameters in line\n@_\n";
                  exit;
               }
               $tmp_string = sprintf "%d %d   %d %d   %d   %e %e %e ;  %s-%s  %s-%s\n",
                                     $atom_i-1, $atom_j-1,  
                                     $atomResNr[$atom_i]-1, $atomResNr[$atom_j]-1,
                                     $funct, $b0, $D, $beta,
                                     $atomName[$atom_i], $atomName[$atom_j],
                                     $atomType[$atom_i], $atomType[$atom_j];
               $bond_info[$mol] .= $tmp_string;
            }
#
#           Harmonic potential
#
            if ( $funct == 1 )
            {
               if ( defined $_[3] and defined $_[4] )
               {
                  $b0 = $_[3];
                  $kb = $_[4];
               }
               else
               {
                  $bond_found = 0;
             
                  if ( $use_opls==1 )
                  {
                     open ITP,"itp_files/ffoplsaabd.itp"
                        or die "Error: cannot open itp file ffoplsaabd.itp !\n";
                  }
                  elsif ( $use_amber03==1 )
                  {
                     open ITP,"itp_files/ffamber03bd.itp"
                        or die "Error: cannot open itp file ffamber03bd.itp !\n";
                  }
                  else
                  {
                     die "Error: cannot find bond parameter for atoms $atom_i $atom_j !\n";
                  }
             
                  while ( <ITP> )
                  {
                     if ( /^\[\s*bondtypes\s*\]/ )
                     {
                        while ( <ITP> )
                        {
#                          skip comments, # lines and blank lines
                           next if ( /^;/ or /^\#/ or /^\s*$/ );
                           next if ( /^\[\s*bondtypes\s*\]/ );
#                          do not read next block
                           last if( /^\[/ );
#                          read bond i, j, funct
                           $_ = (split ";",$_)[0];
                           @_ = (split);
                           $ff_atom_i = $_[0];
                           $ff_atom_j = $_[1];
                           $ff_funct = $_[2];
                           if ( ( $atomType[$atom_i] eq $ff_atom_i 
                                  and $atomType[$atom_j] eq $ff_atom_j 
                                  and $funct eq $ff_funct )
                                or 
                                ( $atomType[$atom_j] eq $ff_atom_i
                                  and $atomType[$atom_i] eq $ff_atom_j
                                  and $funct eq $ff_funct ) )
                           {
                              $b0 = $_[3];
                              $kb = $_[4];
                              $bond_found = 1;
                              last;
                           }
                        }
                        last if $bond_found==1;
                     }
                  }
                  close ITP;
                  die "Error: cannot find bond info for atoms $atom_i $atomType[$atom_i] $atom_j $atomType[$atom_j]!\n"
                     unless $bond_found==1;
               }
               $tmp_string = sprintf "%d %d   %d %d   %d   %e %e  ;  %s-%s  %s-%s\n",
                                     $atom_i-1, $atom_j-1, 
                                     $atomResNr[$atom_i]-1, $atomResNr[$atom_j]-1,
                                     $funct, $b0, $kb,
                                     $atomName[$atom_i], $atomName[$atom_j],
                                     $atomType[$atom_i], $atomType[$atom_j];
               $bond_info[$mol] .= $tmp_string;
            }
         }
      }

#
#     [ paris ] block: scaled 1-4 interaction
#
      if ( /^\[\s*pairs\s*\]/ )
      {
         while ( <TOP> )
         {
#           skip comments, # lines and blank lines
            next if ( /^;/ or /^\#/ or /^\s*$/ );
            next if ( /^\[\s*pairs\s*\]/ );
#           do not read next block
            last if( /^\[/ );
#           read bond i, j, funct
            $_ = (split ";",$_)[0];
            @_ = (split);
            $atom_i = $_[0];
            $atom_j = $_[1];

            #next if $atomResNr[$atom_i]==$atomResNr[$atom_j];
            $nPairs[$mol]++;

            $funct = $_[2];
            if ( $funct == 1 )
            {
               $pair = "";
            }
            elsif ( $funct == 2 )
            {
               shift @_; shift @_; shift @_;
               $pair = sprintf "@_";
            }
            else
            {
               die "Error: unsupported pair function type: $funct !\n";
            }

            $tmp_string = sprintf "%d %d   %d %d   %d   %s  ;  %s-%s  %s-%s\n",
                                  $atom_i-1, $atom_j-1, 
                                  $atomResNr[$atom_i]-1, $atomResNr[$atom_j]-1,
                                  $funct, $pair,
                                  $atomName[$atom_i], $atomName[$atom_j],
                                  $atomType[$atom_i], $atomType[$atom_j];
            $pair_info[$mol] .= $tmp_string;
         }

      }


#     [ angles ] block: angle bending
      if ( /^\[\s*angles\s*\]/ )
      {
         while ( <TOP> )
         {
#           skip comments, # lines and blank lines
            next if ( /^;/ or /^\#/ or /^\s*$/ );
            next if ( /^\[\s*angles\s*\]/ );
#           do not read next block
            last if( /^\[/ );
#           read bond i, j, funct
            $_ = (split ";",$_)[0];
            @_ = (split);
            $atom_i = $_[0];
            $atom_j = $_[1];
            $atom_k = $_[2];

            #next if ( $atomResNr[$atom_i]==$atomResNr[$atom_j]
            #          and $atomResNr[$atom_j]==$atomResNr[$atom_k] );
            $nAngles[$mol]++;

            $funct = $_[3];
            die "Error: unsupported angle function type: $funct !\n"
               unless $funct==1;
            if ( defined $_[4] and defined $_[5] )
            {
               $a0 = $_[4];
               $ka = $_[5];
            }
            else
            {
               $angle_found = 0;

               if ( $use_opls==1 )
               {
                  open ITP,"itp_files/ffoplsaabd.itp"
                     or die "Error: cannot open itp file ffoplsaabd.itp !\n";
               }
               elsif ( $use_amber03==1 )
               {
                  open ITP,"itp_files/ffamber03bd.itp"
                     or die "Error: cannot open itp file ffamber03bd.itp !\n";
               }
               else
               {
                  die "Error: cannot find angle parameter for atoms $atom_i $atom_j $atom_k!\n";
               }

               while ( <ITP> )
               {
                  if ( /^\[\s*angletypes\s*\]/ )
                  {
                     while ( <ITP> )
                     {
#                       skip comments, # lines and blank lines
                        next if ( /^;/ or /^\#/ or /^\s*$/ );
                        next if ( /^\[\s*angletypes\s*\]/ );
#                       do not read next block
                        last if( /^\[/ );
#                       read bond i, j, funct
                        $_ = (split ";",$_)[0];
                        @_ = (split);
                        $ff_atom_i = $_[0];
                        $ff_atom_j = $_[1];
                        $ff_atom_k = $_[2];
                        $ff_funct = $_[3];
                        if ( ( $atomType[$atom_i] eq $ff_atom_i 
                               and $atomType[$atom_j] eq $ff_atom_j 
                               and $atomType[$atom_k] eq $ff_atom_k 
                               and $funct eq $ff_funct )
                             or 
                             ( $atomType[$atom_k] eq $ff_atom_i
                               and $atomType[$atom_j] eq $ff_atom_j
                               and $atomType[$atom_i] eq $ff_atom_k
                               and $funct eq $ff_funct ) )
                        {
                           $a0 = $_[4];
                           $ka = $_[5];
                           $angle_found = 1;
                           last;
                        }
                     }
                     last if $angle_found==1;
                  }
               }
               close ITP;
               die "Error: cannot find angle info for atoms $atom_i $atomType[$atom_i] $atom_j $atomType[$atom_j] $atom_k $atomType[$atom_k] !\n"
                  unless $angle_found==1;
            }
            $tmp_string = sprintf "%d %d %d   %d %d %d   %d   %e %e  ;  %s-%s-%s  %s-%s-%s\n",
                                  $atom_i-1, $atom_j-1, $atom_k-1, 
                                  $atomResNr[$atom_i]-1, $atomResNr[$atom_j]-1, $atomResNr[$atom_k]-1,
                                  $funct, $a0, $ka,
                                  $atomName[$atom_i], $atomName[$atom_j], $atomName[$atom_k],
                                  $atomType[$atom_i], $atomType[$atom_j], $atomType[$atom_k];
            $angle_info[$mol] .= $tmp_string;
         }

      }


#     [ dihedrals ] block: torsion parameters
      if ( /^\[\s*dihedrals\s*\]/ )
      {
         while ( <TOP> )
         {
#           skip comments, # lines and blank lines
            next if ( /^;/ or /^\#/ or /^\s*$/ );
            next if ( /^\[\s*dihedrals\s*\]/ );
#           do not read next block
            last if( /^\[/ );
#           read bond i, j, funct
            $_ = (split ";",$_)[0];
            @_ = (split);
            $atom_i = $_[0];
            $atom_j = $_[1];
            $atom_k = $_[2];
            $atom_l = $_[3];

            #next if ( $atomResNr[$atom_i]==$atomResNr[$atom_j]
            #          and $atomResNr[$atom_j]==$atomResNr[$atom_k] 
            #          and $atomResNr[$atom_k]==$atomResNr[$atom_l] );

            $funct = $_[4];
            die "Error: unsupported dihedral function type: $funct !\n"
               unless ($funct==3 or $funct==1 or $funct==4 or $funct==9);

            if ( defined $_[5] )
            {
               shift @_; shift @_; shift @_; shift @_; shift @_;
               $dih = sprintf "@_";
               $dihedral_found = 1;

#              OPLS proper dihedral
               if ( $funct == 3 )
               {
                  $potential = 0;
                  @coeff = (split " ",$dih);
                  for $ii ( 0 .. (@coeff-1) )
                  {
                     $potential = 1 if $coeff[$ii]!=0.0;
                  }
                  next if $potential==0;
               }
#              OPLS improper dihedral
               elsif ( $funct==1 )
               {
                  if ( $dih =~ /improper/ )
                  {
                     $dih = "180.0 43.93200 2" if ( $dih=~/O_C_X_Y/);
                     $dih = "180.0 43.93200 2" if ( $dih=~/X_NO_ON_NO/);
                     $dih = "180.0 43.93200 2" if ( $dih=~/N2_X_N2_N2/);
                     $dih = "180.0 4.18400 2" if ( $dih=~/Z_N_X_Y/);
                     $dih = "180.0 62.76000 2" if ( $dih=~/Z_CM_X_Y/);
                     $dih = "180.0 4.60240 2" if ( $dih=~/Z_CA_X_Y/);
                  }
                  else
                  {
                     @coeff = (split " ",$dih);
                     next if $coeff[1]==0.0;
                  }
               }
#              Amber03 proper/improper dihedral
               elsif ( $funct==9 or $funct==4 )
               {
                  @coeff = (split " ",$dih);
                  next if $coeff[1]==0.0;
               }
               $nDihedrals[$mol]++;

               $tmp_string = sprintf "%d %d %d %d   %d %d %d %d   %d   %s  ;  %s-%s-%s-%s  %s-%s-%s-%s\n",
                                     $atom_i-1, $atom_j-1, $atom_k-1, $atom_l-1, 
                                     $atomResNr[$atom_i]-1, $atomResNr[$atom_j]-1, 
                                     $atomResNr[$atom_k]-1, $atomResNr[$atom_l]-1,
                                     $funct, $dih, 
                                     $atomName[$atom_i], $atomName[$atom_j], 
                                     $atomName[$atom_k], $atomName[$atom_l],
                                     $atomType[$atom_i], $atomType[$atom_j], 
                                     $atomType[$atom_k], $atomType[$atom_l];
               $dihedral_info[$mol] .= $tmp_string;

            }
            else
            {
               $dihedral_found = 0;

               if ( $use_opls==1 )
               {
                  open ITP,"itp_files/ffoplsaabd.itp"
                     or die "Error: cannot open itp file ffoplsaabd.itp !\n";
               }
               elsif ( $use_amber03==1 )
               {
                  open ITP,"itp_files/ffamber03bd.itp"
                     or die "Error: cannot open itp file ffamber03bd.itp !\n";
               }
               else
               {
                  die "Error: cannot find dihedral parameter for atoms $atom_i $atom_j $atom_k $atom_l!\n";
               }

               while ( <ITP> )
               {
                  if ( /^\[\s*dihedraltypes\s*\]/ )
                  {
                     while ( <ITP> )
                     {
#                       skip comments, # lines and blank lines
                        next if ( /^;/ or /^\#/ or /^\s*$/ );
                        next if ( /^\[\s*dihedraltypes\s*\]/ );
#                       do not read next block
                        last if( /^\[/ );
#                       read bond i, j, funct
                        $_ = (split ";",$_)[0];
                        @_ = (split);
                        $ff_atom_i = $_[0];
                        $ff_atom_j = $_[1];
                        $ff_atom_k = $_[2];
                        $ff_atom_l = $_[3];
                        $ff_funct = $_[4];
                        if ( ( $atomType[$atom_i] eq $ff_atom_i 
                               and $atomType[$atom_j] eq $ff_atom_j 
                               and $atomType[$atom_k] eq $ff_atom_k 
                               and $atomType[$atom_l] eq $ff_atom_l 
                               and $funct eq $ff_funct )
                             or 
                             ( $atomType[$atom_l] eq $ff_atom_i
                               and $atomType[$atom_k] eq $ff_atom_j
                               and $atomType[$atom_j] eq $ff_atom_k
                               and $atomType[$atom_i] eq $ff_atom_l
                               and $funct eq $ff_funct ) )
                        {
                           shift @_; shift @_; shift @_; shift @_; shift @_;
                           $dih = sprintf "@_";
                           $dihedral_found = 1;

                           if ( $funct == 3 )
                           {
                              $potential = 0;
                              @coeff = (split " ",$dih);
                              for $ii ( 0 .. (@coeff-1) )
                              {
                                 $potential = 1 if $coeff[$ii]!=0.0;
                              }
                              next if $potential==0;
                           }
                           elsif ( $funct==1 or $funct==4 or $funct==9 )
                           {
                              @coeff = (split " ",$dih);
                              next if $coeff[1]==0.0;
                           }
                           $nDihedrals[$mol]++;

                           $tmp_string = sprintf "%d %d %d %d   %d %d %d %d   %d   %s  ;  %s-%s-%s-%s  %s-%s-%s-%s\n",
                                                 $atom_i-1, $atom_j-1, $atom_k-1, $atom_l-1, 
                                                 $atomResNr[$atom_i]-1, $atomResNr[$atom_j]-1, 
                                                 $atomResNr[$atom_k]-1, $atomResNr[$atom_l]-1,
                                                 $funct, $dih, 
                                                 $atomName[$atom_i], $atomName[$atom_j], 
                                                 $atomName[$atom_k], $atomName[$atom_l],
                                                 $atomType[$atom_i], $atomType[$atom_j], 
                                                 $atomType[$atom_k], $atomType[$atom_l];
                           $dihedral_info[$mol] .= $tmp_string;

                        }
                     }
                  }
               }
               close ITP;


               if ( $use_amber03==1 and $dihedral_found==0 )
               {
                  open ITP,"itp_files/ffamber03bd.itp"
                     or die "Error: cannot open itp file ffamber03bd.itp !\n";
                  while ( <ITP> )
                  {
                     if ( /^\[\s*dihedraltypes\s*\]/ )
                     {
                        while ( <ITP> )
                        {
#                          skip comments, # lines and blank lines
                           next if ( /^;/ or /^\#/ or /^\s*$/ );
                           next if ( /^\[\s*dihedraltypes\s*\]/ );
#                          do not read next block
                           last if( /^\[/ );
#                          read bond i, j, funct
                           $_ = (split ";",$_)[0];
                           @_ = (split);
                           $ff_atom_i = $_[0];
                           $ff_atom_j = $_[1];
                           $ff_atom_k = $_[2];
                           $ff_atom_l = $_[3];
                           $ff_funct = $_[4];
                           if ( 

                                ( "X" eq $ff_atom_i 
                                  and $atomType[$atom_j] eq $ff_atom_j 
                                  and $atomType[$atom_k] eq $ff_atom_k 
                                  and "X" eq $ff_atom_l 
                                  and $funct eq $ff_funct )
                                or 
                                ( "X" eq $ff_atom_i
                                  and $atomType[$atom_k] eq $ff_atom_j
                                  and $atomType[$atom_j] eq $ff_atom_k
                                  and "X" eq $ff_atom_l
                                  and $funct eq $ff_funct ) 

                                or 

                                ( "X" eq $ff_atom_i 
                                  and $atomType[$atom_j] eq $ff_atom_j 
                                  and $atomType[$atom_k] eq $ff_atom_k 
                                  and $atomType[$atom_l] eq $ff_atom_l 
                                  and $funct eq $ff_funct )
                                or 
                                ( $atomType[$atom_l] eq $ff_atom_i
                                  and $atomType[$atom_k] eq $ff_atom_j
                                  and $atomType[$atom_j] eq $ff_atom_k
                                  and "X" eq $ff_atom_l
                                  and $funct eq $ff_funct ) 

                                or 

                                ( "X" eq $ff_atom_i 
                                  and "X" eq $ff_atom_j 
                                  and $atomType[$atom_k] eq $ff_atom_k 
                                  and $atomType[$atom_l] eq $ff_atom_l 
                                  and $funct eq $ff_funct )
                                or 
                                ( $atomType[$atom_l] eq $ff_atom_i
                                  and $atomType[$atom_k] eq $ff_atom_j
                                  and "X" eq $ff_atom_k
                                  and "X" eq $ff_atom_l
                                  and $funct eq $ff_funct ) 

                              )
                           {
                              shift @_; shift @_; shift @_; shift @_; shift @_;
                              $dih = sprintf "@_";
                              $dihedral_found = 1;

                              if ( $funct == 3 )
                              {
                                 $potential = 0;
                                 @coeff = (split " ",$dih);
                                 for $ii ( 0 .. @coeff )
                                 {
                                    $potential = 1 if $coeff[$ii]!=0.0;
                                 }
                                 next if $potential==0;
                              }
                              elsif ( $funct==1 or $funct==4 or $funct==9 )
                              {
                                 @coeff = (split " ",$dih);
                                 next if $coeff[1]==0.0;
                              }
                              $nDihedrals[$mol]++;

                              $tmp_string = sprintf "%d %d %d %d   %d %d %d %d   %d   %s  ;  %s-%s-%s-%s  %s-%s-%s-%s\n",
                                                    $atom_i-1, $atom_j-1, $atom_k-1, $atom_l-1, 
                                                    $atomResNr[$atom_i]-1, $atomResNr[$atom_j]-1, 
                                                    $atomResNr[$atom_k]-1, $atomResNr[$atom_l]-1,
                                                    $funct, $dih, 
                                                    $atomName[$atom_i], $atomName[$atom_j], 
                                                    $atomName[$atom_k], $atomName[$atom_l],
                                                    $atomType[$atom_i], $atomType[$atom_j], 
                                                    $atomType[$atom_k], $atomType[$atom_l];
                              $dihedral_info[$mol] .= $tmp_string;

                           }
                        }
                     }
                  }
               }
               close ITP;

               die "Error: cannot find dihedral info for atoms $atom_i $atomType[$atom_i] $atom_j $atomType[$atom_j] $atom_k $atomType[$atom_k] $atom_l $atomType[$atom_l] !\n"
                  unless $dihedral_found==1;

            }

         }

      }

   }
   close TOP;

#  add one more molecule type and two more atom types for SPC/E water
   if ( $use_spce==1 ) 
   {
      $mol++;
      $nAtoms = 3;
      @atomCharge = qw/* -0.8476 0.4238 0.4238/;
      @atomMass = qw/* 15.99940 1.00800 1.00800/;
     
#     Oxygen in SPC/E water: 3.16557e-01  6.50194e-01
      $nTypes++;
      $typeName[$nTypes] = 'OW_SPC';
#     add type ID to atom
      $thisType = $nTypes;
      $typeID[1] = $thisType;
      $typeSIG[$typeID[1]] = 0.316557;
      $typeEPS[$typeID[1]] = 0.650194;
     
#     Hydrogen in SPC/E water
      $nTypes++;
      $typeName[$nTypes] = 'HW_SPC';
#     add type ID to atom
      $thisType = $nTypes;
      $typeID[2] = $thisType;
      $typeID[3] = $thisType;
      $typeSIG[$typeID[2]] = 0.0;
      $typeEPS[$typeID[2]] = 0.0;
      $typeSIG[$typeID[3]] = 0.0;
      $typeEPS[$typeID[3]] = 0.0;
     
#     write information of a molecule
#     Note: $typeID[$atom]-1 is printed so that it starts from 0.
#           This is beneficial for C program to read.

      $tmp_string = sprintf "%-15s%-10d%-10d\n", 
                    $molName[$mol], $molNumber[$mol], 1;
      $mol_info .= $tmp_string;

      $tmp_string = sprintf "%-20s%-10d\n", $molName[$mol]."_ATOMS", $nAtoms;
      $atom_info .= $tmp_string;
      for $atom (1..$nAtoms)
      {
         $tmp_string = sprintf "%13.6f%13.6f%15.8f%15.8f%8d%8d\n",
                       $atomCharge[$atom], $atomMass[$atom],
                       $typeSIG[$typeID[$atom]], $typeEPS[$typeID[$atom]],
                       $typeID[$atom]-1, 0;
         $atom_info .= $tmp_string;
      }

#     bonded parameters for SPC/E water model
      $nBonds[$mol] = 2;
      $bond_info[$mol] = "0 1   0 0   1   1.000000e-01 3.450000e+05  ;  OW-HW1  OW-HW\n";
      $bond_info[$mol] .= "0 2   0 0   1   1.000000e-01 3.450000e+05  ;  OW-HW2  OW-HW\n";
      $nPairs[$mol] = 0;
      $pair_info[$mol] = "";
      $nAngles[$mol] = 1;
      $angle_info[$mol] = "1 0 2   0 0 0   1   1.094700e+02 3.830000e+02  ;  HW1-OW-HW2  HW-OW-HW\n";
      $nDihedrals[$mol] = 0;
      $dihedral_info[$mol] = "";

   }

#  check total number of molecules
   die "Error: incorrect number of molecules!\n"
      unless ( $mol==$nMols );

#  write atom info on screen
   printf "Atom_Types   $nTypes\n";
   for $iType (1..$nTypes)
   {
      printf "%s ", $typeName[$iType];
   }
   printf "\n";
   printf "\n";

#  write Molecule and Atom info into param.txt
   printf PAR "$mol_info";
   printf PAR "$atom_info";

#  write Lennard-Jones C6 and C12 coefficients into param.txt
   my $count = 0;
   my $combSIG; my $combEPS;
   printf PAR "%-15s%-10d%-10d\n", "LJ_C6_C12", $nTypes, $nTypes*($nTypes+1)/2;
   for $iType (1..$nTypes)
   {
      for $jType ($iType..$nTypes)
      {
         if ( $comb_rule==2 )
         {
            $combSIG = 0.5*( $typeSIG[$iType] + $typeSIG[$jType] );
            $combEPS = sqrt( $typeEPS[$iType] * $typeEPS[$jType] );
         }
         elsif ( $comb_rule==3 )
         {
            $combSIG = sqrt( $typeSIG[$iType] * $typeSIG[$jType] );
            $combEPS = sqrt( $typeEPS[$iType] * $typeEPS[$jType] );
         }
         else
         {
            die "Error: unrecognized comb_rule: $comb_rule !\n";
         }
         $count++;
#        write parameters
#        Note: iType-1, jType-1 and count-1 are printed so that
#              they all start from 0.
         printf PAR "%20.8e%20.8e%8d%8d%8d\n", 
            4.0 * $combEPS * $combSIG**6, 
            4.0 * $combEPS * $combSIG**12,
            $iType-1, $jType-1, $count-1;
#
#           Note: $count-1 = $nTypes*($nTypes+1)/2 
#                          - ($nTypes-($iType-1))*($nTypes-($iType-1)+1)/2 
#                          + ($jType-1)+1-($iType-1)-1
#                          = ($nTypes*2-1-($iType-1))*($iType-1)/2+($jType-1)
#
      }
   }

#  write bonded potentials into param.txt
   for $mol ( 1 .. $nMols )
   {
      $tmp_string = sprintf "%-20s%-10d\n", $molName[$mol]."_BONDS", $nBonds[$mol];
      $bond_info[$mol] = $tmp_string.$bond_info[$mol];
      $tmp_string = sprintf "%-20s%-10d\n", $molName[$mol]."_PAIRS", $nPairs[$mol];
      $pair_info[$mol] = $tmp_string.$pair_info[$mol];
      $tmp_string = sprintf "%-20s%-10d\n", $molName[$mol]."_ANGLES", $nAngles[$mol];
      $angle_info[$mol] = $tmp_string.$angle_info[$mol];
      $tmp_string = sprintf "%-20s%-10d\n", $molName[$mol]."_DIHEDRALS", $nDihedrals[$mol];
      $dihedral_info[$mol] = $tmp_string.$dihedral_info[$mol];

      printf PAR "$bond_info[$mol]";
      printf PAR "$pair_info[$mol]";
      printf PAR "$angle_info[$mol]";
      printf PAR "$dihedral_info[$mol]";

#     write bonded info on screen
      printf "%-20s%-10d\n", $molName[$mol]."_BONDS", $nBonds[$mol];
      printf "%-20s%-10d\n", $molName[$mol]."_PAIRS", $nPairs[$mol];
      printf "%-20s%-10d\n", $molName[$mol]."_ANGLES", $nAngles[$mol];
      printf "%-20s%-10d\n", $molName[$mol]."_DIHEDRALS", $nDihedrals[$mol];
      printf "\n";

   }

   printf PAR "END\n";

   printf "Data written to file param.txt\n";
   printf "\n";

#  end of program

