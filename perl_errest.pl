=comment
 This file is part of the SurfTensCalc program.                                      
 Copyright (C) 2012 Xin Li
                                                                           
 Filename:  perl_errest.pl
 License:   BSD 2-Clause License

 This software is provided by the copyright holders and contributors "as is"
 and any express or implied warranties, including, but not limited to, the
 implied warranties of merchantability and fitness for a particular purpose are
 disclaimed. In no event shall the copyright holder or contributors be liable
 for any direct, indirect, incidental, special, exemplary, or consequential
 damages (including, but not limited to, procurement of substitute goods or
 services; loss of use, data, or profits; or business interruption) however
 caused and on any theory of liability, whether in contract, strict liability,
 or tort (including negligence or otherwise) arising in any way out of the use
 of this software, even if advised of the possibility of such damage.
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


