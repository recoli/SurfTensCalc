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

open LOG,">log" or die;

unless (@ARGV==4) {
  printf LOG "Incorrect number of arguments!\n";
  exit(1);
}
$nFrames = $ARGV[0];
$nStart  = $ARGV[1];
$filenameXTC = $ARGV[2];
$filenameTOP = $ARGV[3];

print LOG "\n";
print LOG "***********************************************\n";
print LOG "*         Surface Tension Program             *\n";
print LOG "*          Written in C and Perl              *\n";
print LOG "*                by Xin Li                    *\n";
print LOG "*          TheoChem, KTH, Sweden              *\n";
print LOG "***********************************************\n";
print LOG "\n";

printf LOG "total number of frames: %-10d\n", $nFrames;
printf LOG "starting from frame:    %-10d\n", $nStart;
print LOG "\n";

printf LOG "trajectory file name:   %-20s\n", $filenameXTC;
printf LOG "topology file name:     %-20s\n", $filenameTOP;
printf LOG "force field file name:  %-20s\n", "ffoplsaanb.itp";
print LOG "\n";

unless (-e "$filenameXTC") {
  printf LOG "Cannot find file $filenameXTC!\n";
  exit(1);
}

unless (-e "$filenameTOP") {
  printf LOG "Cannot find file $filenameTOP!\n";
  exit(1);
}

unless (-e "ffoplsaanb.itp") {
  printf LOG "Cannot find file ffoplsaanb.itp!\n";
  exit(1);
}

#******************************************
# Run perl_top 
#******************************************

print LOG "Parsing topology files and writing to param.txt ... ";
system "perl perl_top.pl $nFrames $nStart $filenameXTC $filenameTOP";
unless (-e "param.txt") {
  printf LOG "Failed to generate param.txt from perl_top.pl !\n";
  exit(1);
}
print LOG "Done\n\n";

#************************************************
# Compile and run the program calc_pres_dens.c
#************************************************

print LOG "Calculating pressure profile and work of formation (most time-consuming) ... ";
system "./calc_pres_dens.x";
unless (-e "pressure.dat") {
  printf LOG "Failed to calculate pressure.dat by using omp_calc_pres_dens.x !\n";
  exit(1);
}
unless (-e "work.dat") { 
  printf LOG "Failed to calculate work.dat by using omp_calc_pres_dens.x !\n";
  exit(1);
}
print LOG "Done\n\n";

#************************************************
# Read data.dat and calculate surface tension
#************************************************

print LOG "Reading output file perssure.dat ... ";

## constants ##

$pi     = 3.14159265 ;
$kT     = 1.380651 * 298.0 * 6.022142 * 0.001 ; # at 298.0 K

$maxbin = 300 ;
$dr     = 0.03;

## read data ##

# r = radius, d = density, p = pressure_U
$bin = 0;
open DAT,"pressure.dat" or die;
while(<DAT>){
  if(/^\#/ or /^\@/){ next; }
  $r[$bin] = (split)[0];
  $d[$bin] = (split)[1];
  $p[$bin] = (split)[2];
  $bin++;
}
close DAT;

print LOG "Done\n\n";

#******************************************************
# fitting density profile to tanh function
#******************************************************

print LOG "Fitting density ... ";

system "./fit_dens.x";
unless (-e "fit_result") {
  printf LOG "Failed to get fit_result from compile fit_dens.x !\n";
  exit(1);
}

print LOG "Done\n\n";

#***************************************
# Calculate surface tension 
#***************************************

print LOG "Calculating surface tension ... ";

open RESULT,"fit_result" or die;
while(<RESULT>){
  if(/Liquid Density/){ $fit_dens = (split)[3]; }
  if(/Re/){ $fit_re = (split)[2]; }
}
close RESULT;

# calculate work of formation
open ST,">surftens.dat" or die;
$work = 0.0;
for $bin (1..$maxbin) {
  $work += ( ($kT*$d[$bin]+$p[$bin])*$r[$bin]**2 ) * $dr;
}
$work *=  $pi * 2.0 ;  # kJ/mol
$surftens =  3.0 * $work / ( 4.0 * $pi * $fit_re ** 2 ) ; # kJ/mol/nm^2
$surftens /=  0.60221415 ;  # mN/m
$work /=  60.221415  ;  # 10^-19 J

printf ST "Density of the cluster     = %13.6F nm^-3\n",     $fit_dens;
printf ST "Equimolar dividing surface = %13.6F nm\n",        $fit_re;
printf ST "Work of formation          = %13.6F e-19 J\n",    $work;
printf ST "Effective surface tension  = %13.6F mN m^-1\n",   $surftens;
close ST;

print LOG "Done\n\n";
print LOG "All Done!\n\n";
close LOG;

#******************
# end of program 
#******************

