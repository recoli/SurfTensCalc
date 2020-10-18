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

#  print title
   printf "\n";
   printf "***********************************************\n";
   printf "*      Perl script for gmx top file           *\n";
   printf "*                by Xin Li                    *\n";
   printf "*          TheoChem, KTH, Sweden              *\n";
   printf "*           Version 2012-01-07                *\n";
   printf "***********************************************\n";

   printf "\n";
   printf "Usage: perl perl_top.pl \\\n";
   printf "       -top [name_of_gmx_top_file] \\\n";
   printf "       -xtc [name_of_gmx_xtc_file] \n";
   printf "\n";

#  options
#  Note: change array to hash
   my %options = @ARGV;

   die "Error: undefined -top option!\n" unless (defined $options{"-top"});
   my $gmx_top = $options{"-top"};
   die "Error: cannot find file $gmx_top !\n" unless (-e $gmx_top);

   die "Error: undefined -xtc option!\n" unless (defined $options{"-xtc"});
   my $gmx_xtc = $options{"-xtc"};
   die "Error: cannot find file $gmx_xtc !\n" unless (-e $gmx_xtc);

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
#           do not read next block
            last if( /^\[/ );
#           print to file
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
         unless (-e "ffoplsaanb.itp") ;
      printf "Will use the OPLS force field. (ffoplsaanb.itp)\n";
   }
   if ( $use_amber03 == 1 )
   {
      die "Error: cannot find file ffamber03nb.itp!\n"
         unless (-e "ffamber03nb.itp") ;
      printf "Will use the AMBER03 force field. (ffamber03nb.itp)\n";
   }
   if ( $use_spce == 1 )
   {
      printf "Will use the SPC/E water model.\n";
   }

   printf "The combination rule for non-bonded LJ is $comb_rule.\n";
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
   printf PAR "%-15s%-20s\n", "xtc_file", $gmx_xtc;

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

#  print molecule info on screen
   printf "Found $nMols molecule types\n";
   for $mol (1..$nMols)
   {
      printf "%s\n", $molName[$mol];
   }
   printf "\n";

#  read and write atomic information
   open TOP,"$gmx_top" 
      or die "Error: cannot open gmx top file $gmx_top !\n";
#  initialize number of types and number of molecules
   my $iType; my $jType; my $nTypes; 
   my $isNewType; my $thisType;
   my $atom; my $nAtoms;

   my @atomResNr; my @atomResName; 
   my @atomName; my @atomType; my @atomCharge; my @atomMass;

   my @typeName; my @typeSIG; my @typeEPS;

   my $typeFound;

   $mol = 0;
   $nTypes = 0;
   while ( <TOP> )
   {
      if ( /^\[\s*atoms\s*\]/ )
      {
         $mol++;
         while ( <TOP> )
         {
#           skip comments or blank lines
            next if ( /^;/ or /^\s*$/ );
#           do not read next block
            last if( /^\[/ );

#           read atom ID, resnr, resname, name, charge, mass
            $atom = (split)[0];
            $nAtoms = $atom;
            $atomResNr[$atom]   = (split)[2];
            $atomResName[$atom] = (split)[3];
            $atomName[$atom]    = (split)[4];
            $atomCharge[$atom]  = (split)[6];
            $atomMass[$atom]    = (split)[7];

#           check if this type is a new type 
            $isNewType = 1;
            for $iType (1..$nTypes) 
            {
#              NOTE: use "eq" to compare two strings (not "==")!
               if ( $typeName[$iType] eq (split)[1] ) 
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
               $typeName[$nTypes] = (split)[1];

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
                  open ITP,"ffoplsaanb.itp"
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
                  open ITP,"ffamber03nb.itp"
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
            $atomType[$atom] = $thisType;

         }
#        printf information of a molecule
#        Note: $atomType[$atom]-1 is printed so that it starts from 0.
#              This is beneficial for C program to read.

         $tmp_string = sprintf "%-15s%-10d%-10d\n", 
                       $molName[$mol], $molNumber[$mol], $atomResNr[$nAtoms];
         $mol_info .= $tmp_string;

         $tmp_string = sprintf "%-20s%-10d\n", $molName[$mol]."_ATOMS", $nAtoms;
         $atom_info .= $tmp_string;
         for $atom (1..$nAtoms)
         {
            $tmp_string = sprintf "%13.6f%13.6f%15.8f%15.8f%6d%6d\n",
                          $atomCharge[$atom], $atomMass[$atom],
                          $typeSIG[$atomType[$atom]], $typeEPS[$atomType[$atom]],
                          $atomType[$atom]-1, $atomResNr[$atom]-1;
            $atom_info .= $tmp_string;
         }
      }
   }
   close TOP;

#  add one more molecule type and two more atom types for spc/e water
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
      $atomType[1] = $thisType;
      $typeSIG[$atomType[1]] = 0.316557;
      $typeEPS[$atomType[1]] = 0.650194;
     
#     Hydrogen in SPC/E water
      $nTypes++;
      $typeName[$nTypes] = 'HW_SPC';
#     add type ID to atom
      $thisType = $nTypes;
      $atomType[2] = $thisType;
      $atomType[3] = $thisType;
      $typeSIG[$atomType[2]] = 0.0;
      $typeEPS[$atomType[2]] = 0.0;
      $typeSIG[$atomType[3]] = 0.0;
      $typeEPS[$atomType[3]] = 0.0;
     
#     printf information of a molecule
#     Note: $atomType[$atom]-1 is printed so that it starts from 0.
#           This is beneficial for C program to read.

      $tmp_string = sprintf "%-15s%-10d%-10d\n", 
                    $molName[$mol], $molNumber[$mol], 1;
      $mol_info .= $tmp_string;

      $tmp_string = sprintf "%-20s%-10d\n", $molName[$mol]."_ATOMS", $nAtoms;
      $atom_info .= $tmp_string;
      for $atom (1..$nAtoms)
      {
         $tmp_string = sprintf "%13.6f%13.6f%15.8f%15.8f%6d%6d\n",
                       $atomCharge[$atom], $atomMass[$atom],
                       $typeSIG[$atomType[$atom]], $typeEPS[$atomType[$atom]],
                       $atomType[$atom]-1, 0;
         $atom_info .= $tmp_string;
      }
   }

   print PAR "$mol_info";
   print PAR "$atom_info";

#  check total number of molecules
   die "Error: incorrect number of molecules!\n"
      unless ( $mol==$nMols );

#  print atom info on screen
   printf "Found $nTypes atom types\n";
   for $iType (1..$nTypes)
   {
      printf "%s\n", $typeName[$iType];
   }
   printf "\n";

#  print Lennard-Jones C6 and C12 coefficients
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
#        print data
#        Note: iType-1, jType-1 and count-1 are printed so that
#              they all start from 0.
         printf PAR "%20.8e%20.8e%6d%6d%6d\n", 
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

#  end of program

   printf PAR "END\n";
   close TOP;
   printf "Data written to file param.txt\n";

   printf "\n";
   printf "<NOTE> Please center the trajectory before surface tension calculations:\n";
   printf "echo 0 0|trjconv -f traj.xtc -s topol.tpr -o new-traj.xtc -ndec 4 -pbc mol -center\n";
   printf "mv new-traj.xtc traj.xtc\n";
   printf "\n";

