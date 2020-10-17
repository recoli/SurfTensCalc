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

    die unless @ARGV==4;
    $nFrames = $ARGV[0];
    $nStart  = $ARGV[1];
    $filenameXTC = $ARGV[2];
    $filenameTOP = $ARGV[3];

#   open file for writing
    open PAR,">param.txt" or die;

    printf PAR "%-10s%10d\n", "nFrames", $nFrames;
    printf PAR "%-10s%10d\n", "nStart", $nStart;
    printf PAR "%-10s%10s\n", "xtc", $filenameXTC;

#   read and write molecular information
    open TOP,"$filenameTOP" or die;
    $mol = 0;
    while(<TOP>){
      if(/^\[\s*molecules\s*\]/){
        while(<TOP>){
          if(/^;/ or /^\s*$/){
            next;
          }else{
            $mol++;
            $molName[$mol] = (split)[0];
            $molNumber[$mol] = (split)[1];
          }
        }
      }
    }
    $nMols = $mol;
    printf PAR "MOLECULES%10d\n", $nMols;
    for $mol (1..$nMols) {
      printf PAR "%5s%10d\n", $molName[$mol], $molNumber[$mol];
    }
    close TOP;

#   read and write atomic information
    open TOP,"$filenameTOP" or die;
#   initialize number of types and number of molecules
    $nTypes = 0;
    $nMols = 0;
    while(<TOP>){
      if(/^\[\s*atoms\s*\]/){
        $nMols++;
        while(<TOP>){
#         skip comments or blank lines
          if(/^;/ or /^\s*$/){ 
            next; 
#         do not read next block
          }elsif(/^\[/){ 
            last; 
#         read this block
          }else{
#           read atom ID, resnr, resname, name, charge, mass
            $atom = (split)[0];
            $nAtoms = $atom;
            $atomResNr[$atom]   = (split)[2];
            $atomResName[$atom] = (split)[3];
            $atomName[$atom]    = (split)[4];
            $atomCharge[$atom]  = (split)[6];
            $atomMass[$atom]    = (split)[7];
#           check if this type is a new type 
            $newType = 1;
            for $type (1..$nTypes) {
#             NOTE: To compare two strings you must use "eq" instead of "="!
              if ( $typeName[$type] eq (split)[1] ) { 
                $newType = 0;
                $theType = $type;
                last;
              }
            }
#           add a new type if this is new
            if ( $newType == 1 ) { 
              $nTypes++;
              $typeName[$nTypes] = (split)[1];
              open ITP,"ffoplsaanb.itp" or die;
              while (<ITP>) {
                if (/^\s+$typeName[$nTypes]\s+/) {
                  $typeSigma[$nTypes] = (split)[6];
                  $typeEpsilon[$nTypes] = (split)[7];
                }
              }
              close ITP;
              $theType = $nTypes;
            }
#           add type ID to atom
            $atomType[$atom] = $theType;
          }
        }
#       print information of a molecule
        printf PAR "%-10s%10d\n", $molName[$nMols], $nAtoms;
        for $atom (1..$nAtoms){
          printf PAR "%13.6f%13.6f%15.8f%15.8f\n",
                 $atomCharge[$atom], $atomMass[$atom],
                 $typeSigma[$atomType[$atom]], $typeEpsilon[$atomType[$atom]];
        }
      }
    }

#   check if spce.itp is included
    $spce = 0;
    open TOP,"$filenameTOP" or die;
    while(<TOP>){
      if(/^\#include\s+\".*?spce\.itp\"/){ $spce = 1; }
    }
    close TOP;

#   add one more molecule type and two more atom types for spc/e water */
    if ( $spce==1 ) {
      $nMols++;
      $nAtoms = 3;
      @atomCharge = qw/* -0.8476 0.4238 0.4238/;
      @atomMass = qw/* 15.99940 1.00800 1.00800/;
     
      $nTypes++;
      $typeName[$nTypes] = 'opls_116';
      open ITP,"ffoplsaanb.itp" or die;
      while (<ITP>) {
        if (/$typeName[$nTypes]/) {
          $typeSigma[$nTypes] = (split)[6];
          $typeEpsilon[$nTypes] = (split)[7];
        }
      }
      close ITP;
      $theType = $nTypes;
#     add type ID to atom
      $atomType[1] = $theType;
     
      $nTypes++;
      $typeName[$nTypes] = 'opls_117';
      open ITP,"ffoplsaanb.itp" or die;
      while (<ITP>) {
        if (/$typeName[$nTypes]/) {
          $typeSigma[$nTypes] = (split)[6];
          $typeEpsilon[$nTypes] = (split)[7];
        }
      }
      close ITP;
      $theType = $nTypes;
#     add type ID to atom
      $atomType[2] = $theType;
      $atomType[3] = $theType;
     
      printf PAR "%-10s%10d\n", $molName[$nMols], $nAtoms;
      for $atom (1..$nAtoms){
        printf PAR "%13.6f%13.6f%15.8f%15.8f\n",
               $atomCharge[$atom], $atomMass[$atom],
               $typeSigma[$atomType[$atom]], $typeEpsilon[$atomType[$atom]];
      }
    }

#   end of program

    printf PAR "END\n";
    close TOP;
    
