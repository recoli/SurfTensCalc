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

   my (@ave,@std,@sem);
   my ($st_cut,$st_full);

   unless (-e "work.dat") {
     print "Cannot find file work.dat!\n";
     exit(1);
   }
 
   system "cp work.dat work.xvg";
   system "g_analyze -f work.xvg > work.SEM.dat";
   open FL,"work.SEM.dat" or die;
   while(<FL>)
   {
      if(/SS(\d+)/)
      { 
         @_ = (split);
         $ave[$1] = $_[1]; 
         $std[$1] = $_[2]; 
         $sem[$1] = $_[3]; 
      }
   }
   close FL;
 
   open DAT,"surftens.dat" or die "Cannot open file surftens.dat !\n";
   while(<DAT>)
   {
      if(/Surface tension  gamma/)
      { 
         @_=(split); 
         $st_cut=$_[4]; 
         $st_full=$_[5];
      }
   }
   close DAT;

#         "   Surface tension  gamma  =      xx.xx      xx.xx  ( mJ/m^2 ) "
   printf "   Standard error of mean  = %10.2f %10.2f  ( mJ/m^2 )\n", 
              $st_cut*$sem[2]/$ave[2], $st_full*$sem[1]/$ave[1];
   printf "   \n", 


