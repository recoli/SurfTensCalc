/*
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
 */

/********************************************************************
 *    C program to analyze data from sphere_calc_pres_dens.x        *
 *    by Xin Li, TheoChem & Biology, KTH, Stockholm, Sweden         *
 *                      Version 2012-03-16                          *
 ********************************************************************/

   #include <math.h>
   #include <stdio.h>
   #include <unistd.h>
   #include <stdlib.h>
   #include <string.h>

   #include "elg-fit.h"
   #include "st-calc.h"

/********************************************************************
 *    Main program                                                  *
 ********************************************************************/

   int main( int argc, char *argv[] )
   {
/*
 *    There should be 2 argument: numprocs, boxSize
 */
      if ( argc != 3 )
      {
         printf("%s\n", "Error: Incorrect number of arguments (should be 2)!");
         printf("%s\n", "1st argument: numprocs");
         printf("%s\n", "2nd argument: approx. box size (in nm)");
         exit(1);
      }
/*
 *    number of processors (data files)
 */
      int iproc, numprocs;
      numprocs = atoi(argv[1]);

      int nFrames;

      double boxSize;
      boxSize = atof(argv[2]); /* nm */

      Vec_elg  elg_param;
      Vec_dens_st  dens_st;
      Vec_pres_st  pres_st_full, pres_st_cut;

/*
 *    Define variables
 */

      double dl,dg,re,r0,ksi;
      double pl,pg,rs;

/*
 *    Constants    
 *
 *    Source: 2010 CODATA
 *    nA (Avogadro constant):        6.02214129 * 10^23 mol^-1
 *
 *    PI to 20 decimal places:       3.14159265358979323846
 */
      const double nA   =   6.02214129 ;  /* 10^23 mol^-1 */
      const double pi   =   3.14159265358979323846 ;

      int bin;
      const double dR = 0.05; /*nm*/
      int maxBin ;
      maxBin = (int)( ( boxSize*0.5 - 1.0 ) / dR ) ;

      char filename[20], line[256];

      char text_work[256], text_pres[256], text_pT[265], text_dens[256];
      char text_start[256], text_error[256], text_end[256];
      FILE *file_data, *file_work, *file_pres, *file_dens, *file_xvg;
      int read_text, read_data;

      int error_found, normal_end;


      double r ;
      double *avR ;

      double pN, pNc, pN_K, pN_U_LJ, pN_U_LJc, pN_U_QQ, pN_U_14, pN_U_BD1, pN_U_BD2 ;
      double *avPN, *avPNc, *avPN_K;
      double *avPN_U_LJ, *avPN_U_LJc, *avPN_U_QQ;
      double *avPN_U_14, *avPN_U_BD1, *avPN_U_BD2;

      double pT, pTc, pT_K, pT_U_LJ, pT_U_LJc, pT_U_QQ, pT_U_14, pT_U_BD1, pT_U_BD2 ;
      double *avPT, *avPTc, *avPT_K;
      double *avPT_U_LJ, *avPT_U_LJc, *avPT_U_QQ;
      double *avPT_U_14, *avPT_U_BD1, *avPT_U_BD2;

      int mol, dens, molTypes;
      double *avDens, **avResDens;
      char *tmp;

      char command[256];
      FILE *file_fit_d, *file_fit_p;

      double fit_dl, fit_dg, fit_re;

      double work_cut, eff_st_cut;
      double fit_pl_cut, fit_pg_cut, fit_rs_cut, fit_st_cut, fit_work_cut, deltaP_cut;

      double work_full, eff_st_full;
      double fit_pl_full, fit_pg_full, fit_rs_full, fit_st_full, fit_work_full, deltaP_full;

      double integr, work_N_full;
      double P_beta, P_0;
      double *avDPN;

/*
 *    set initial values and allocate arrays
 */
      molTypes = 5000;

      avR = malloc( sizeof(double) * maxBin );

      avPN = malloc( sizeof(double) * maxBin );
      avPNc = malloc( sizeof(double) * maxBin );
      avPN_K = malloc( sizeof(double) * maxBin );
      avPN_U_LJ = malloc( sizeof(double) * maxBin );
      avPN_U_LJc = malloc( sizeof(double) * maxBin );
      avPN_U_QQ = malloc( sizeof(double) * maxBin );
      avPN_U_14 = malloc( sizeof(double) * maxBin );
      avPN_U_BD1 = malloc( sizeof(double) * maxBin );
      avPN_U_BD2 = malloc( sizeof(double) * maxBin );

      avPT = malloc( sizeof(double) * maxBin );
      avPTc = malloc( sizeof(double) * maxBin );
      avPT_K = malloc( sizeof(double) * maxBin );
      avPT_U_LJ = malloc( sizeof(double) * maxBin );
      avPT_U_LJc = malloc( sizeof(double) * maxBin );
      avPT_U_QQ = malloc( sizeof(double) * maxBin );
      avPT_U_14 = malloc( sizeof(double) * maxBin );
      avPT_U_BD1 = malloc( sizeof(double) * maxBin );
      avPT_U_BD2 = malloc( sizeof(double) * maxBin );

      avDens = malloc( sizeof(double) * maxBin );
      avResDens = malloc(maxBin * sizeof(double *));
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         avResDens[bin] = malloc(molTypes * sizeof(double));
      }

      avDPN = malloc( sizeof( double ) * maxBin );

/*
 *    initialize arrays
 */
      for ( bin=0; bin<maxBin; bin++ )
      {
         avR[bin] = ( (double)(bin) + 0.5 ) * dR;

         avPN[bin] = 0.0;
         avPNc[bin] = 0.0;
         avPN_K[bin] = 0.0;
         avPN_U_LJ[bin] = 0.0;
         avPN_U_LJc[bin] = 0.0;
         avPN_U_QQ[bin] = 0.0;
         avPN_U_14[bin] = 0.0;
         avPN_U_BD1[bin] = 0.0;
         avPN_U_BD2[bin] = 0.0;

         avPT[bin] = 0.0;
         avPTc[bin] = 0.0;
         avPT_K[bin] = 0.0;
         avPT_U_LJ[bin] = 0.0;
         avPT_U_LJc[bin] = 0.0;
         avPT_U_QQ[bin] = 0.0;
         avPT_U_14[bin] = 0.0;
         avPT_U_BD1[bin] = 0.0;
         avPT_U_BD2[bin] = 0.0;

         avDens[bin] = 0.0;
         for ( mol=0; mol<molTypes; mol++ )
         {
            avResDens[bin][mol] = 0.0;
         }
      }
/*
 *    read data files
 */

      printf( "<> Collecting data from each processor...\n" );

      sprintf ( text_work, 
                "#Frame  W_N_full  W_N_cut  W_T_full  W_T_cut\n" );
      sprintf ( text_pres, 
                "#Radius  pN  pNc  pN_K  pN_U_LJ  pN_U_LJc  pN_U_QQ  pN_U_14  pN_U_BD1  pN_U_BD2\n" );
      sprintf ( text_pT,  
                "#Radius  pT  pTc  pT_K  pT_U_LJ  pT_U_LJc  pT_U_QQ  pT_U_14  pT_U_BD1  pT_U_BD2\n" ); 
      sprintf ( text_dens, 
                "#Radius  Total_Density  Residue_Densities (g/cm^3)\n" );

      sprintf ( text_error, "Error:" );
      sprintf ( text_end, ">>> Calculation ended normally\n" );

      file_work = fopen( "work.dat", "w" ) ;

      nFrames = 0;

      for ( iproc=0; iproc<numprocs; iproc++ )
      {

         printf( "   Processor %d ... ", iproc );
/*
 *       gunzip log file
 */     
         sprintf ( filename, "data-%d.log.gz", iproc );

         if( access( filename, F_OK ) == -1 ) {
            printf( "\n<> Error: File %s doesn't exist!\n", filename ) ;
            exit(1);
         }
         sprintf ( command, "gunzip -f %s", filename );
         system ( command );

         sprintf ( filename, "data-%d.log", iproc );
/*
 *       check log file
 */
         file_data = fopen( filename, "r" ) ;
         if ( NULL==file_data )
         {
            printf( "Error: Cannot open file %s !\n", filename ) ;
            exit(1);
         }
         normal_end = 0;
         error_found = 0;
         while( fgets(line,sizeof(line),file_data) != NULL )
         {
            if ( strcmp(line,text_end) == 0 )
            {
               normal_end = 1;
            }
            sscanf( line, "%s", text_start );
            if ( strcmp(text_start,text_error) == 0 )
            {
               error_found = 1;
               printf("\n<> Warning: Error message found in file %s\n", filename);
               printf("   %s", line);
               continue;
            }
         }
         fclose( file_data );
         if ( normal_end == 0 || error_found == 1 )
         {
            if ( error_found == 0 ) { printf("\n"); }
            if ( normal_end == 0 )
            {
               printf("<> Warning: Ending message not found in file %s\n", filename);
            }
            printf("   Skipping file %s\n", filename);
            continue;
         }

/*
 *       read log file
 */     
         file_data = fopen( filename, "r" ) ;
         if ( NULL==file_data ) 
         {
            printf( "Error: Cannot open file %s !\n", filename ) ;
            exit(1);
         }
       
         read_text = 0;
         read_data = 0;
         while( fgets(line,sizeof(line),file_data) != NULL )
         {
            if ( strcmp(line,text_work) == 0 ) 
            {
               read_text = 1;
               read_data = 0;
               fgets( line, sizeof(line), file_data );
               fprintf( file_work, "%s", line );
               nFrames++;
            }
            else if ( strcmp(line,text_pres) == 0 ) 
            {
               read_text = 2;
               read_data = 0;
               for ( bin=0; bin<maxBin; bin++ )
               {
                  fgets( line, sizeof(line), file_data );
                  sscanf( line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
                          &r, &pN, &pNc, &pN_K, 
                          &pN_U_LJ, &pN_U_LJc, &pN_U_QQ, 
                          &pN_U_14, &pN_U_BD1, &pN_U_BD2 );
                  avPN[bin] += pN;
                  avPNc[bin] += pNc;
                  avPN_K[bin] += pN_K;
                  avPN_U_LJ[bin]  += pN_U_LJ;
                  avPN_U_LJc[bin] += pN_U_LJc;
                  avPN_U_QQ[bin]  += pN_U_QQ;
                  avPN_U_14[bin]  += pN_U_14;
                  avPN_U_BD1[bin] += pN_U_BD1;
                  avPN_U_BD2[bin]  += pN_U_BD2;
               }
            }
            else if ( strcmp(line,text_pT) == 0 ) 
            {
               read_text = 3;
               read_data = 0;
               for ( bin=0; bin<maxBin; bin++ )
               {
                  fgets( line, sizeof(line), file_data );
                  sscanf( line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
                          &r, &pT, &pTc, &pT_K, 
                          &pT_U_LJ, &pT_U_LJc, &pT_U_QQ, 
                          &pT_U_14, &pT_U_BD1, &pT_U_BD2 );
                  avPT[bin] += pT;
                  avPTc[bin] += pTc;
                  avPT_K[bin] += pT_K;
                  avPT_U_LJ[bin]  += pT_U_LJ;
                  avPT_U_LJc[bin] += pT_U_LJc;
                  avPT_U_QQ[bin]  += pT_U_QQ;
                  avPT_U_14[bin]  += pT_U_14;
                  avPT_U_BD1[bin] += pT_U_BD1;
                  avPT_U_BD2[bin]  += pT_U_BD2;
               }
            }
            else if ( strcmp(line,text_dens) == 0 ) 
            {
               read_text = 4;
               read_data = 0;
               for ( bin=0; bin<maxBin; bin++ )
               {
                  fgets( line, sizeof(line), file_data );
/*
 *                read the 1st column: radius
 */
                  tmp = strtok( line, " " );
/*
 *                read the 2nd column: total density
 */
                  tmp = strtok( NULL, " " );
                  if ( NULL==tmp )
                  {
                     break;
                  }
                  avDens[bin] += atof(tmp);
/*
 *                read the rest columns: individual densities
 */
                  mol = 0;
                  while( 1 ) 
                  {
                     tmp = strtok( NULL, " " );
                     if ( NULL==tmp ) 
                     {
                        break;
                     }
                     avResDens[bin][mol] += atof(tmp);
                     mol ++;
                  } 
                  molTypes = mol;

               }
            }
            else
            {
               read_data = 1;
            }

/*
            if ( read_text==1 && read_data==1 ) 
            {
               fprintf( file_work, "%s", line );
            }
*/

         }

         fclose( file_data );
/*
 *       gzip log file
 */     
         sprintf ( command, "gzip -f %s", filename );
         system ( command );

         printf( "Done\n" );

      }

      if ( nFrames == 0 )
      {
         printf("<> Error: No frame is read!\n");
         exit(1);
      }

      printf( "<> Number of Frames: %d\n", nFrames );
      printf( "   Printing work of formation to work.dat\n" );

/*
 *    write pressure.dat
 */
      printf( "   Printing pressure_N profile to pressure.dat\n" );

      file_pres = fopen( "pressure_all.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         avPN[bin]  /= nFrames;
         avPNc[bin] /= nFrames;
         avPN_K[bin] /= nFrames;
         avPN_U_LJ[bin]  /= nFrames;
         avPN_U_LJc[bin] /= nFrames;
         avPN_U_QQ[bin]  /= nFrames;
         avPN_U_14[bin]  /= nFrames;
         avPN_U_BD1[bin] /= nFrames;
         avPN_U_BD2[bin] /= nFrames;
         fprintf ( file_pres, "%.3f  %.6f %.6f %.6f  %.6f %.6f %.6f  %.6f %.6f %.6f\n", 
                   avR[bin], avPN[bin], avPNc[bin], avPN_K[bin], 
                   avPN_U_LJ[bin], avPN_U_LJc[bin], avPN_U_QQ[bin],
                   avPN_U_14[bin], avPN_U_BD1[bin], avPN_U_BD2[bin] );
      }
      fclose( file_pres );

      file_pres = fopen( "pressure_cut.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_pres, "%.3f  %.6f %.6f %.6f\n",
                   avR[bin], avPNc[bin], avPN_K[bin], avPNc[bin]-avPN_K[bin] );
      }
      fclose( file_pres );

      file_pres = fopen( "pressure_full.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_pres, "%.3f  %.6f %.6f %.6f\n",
                   avR[bin], avPN[bin], avPN_K[bin], avPN[bin]-avPN_K[bin] );
      }
      fclose( file_pres );
/*
 *    write pT.dat
 */
      printf( "   Printing pressure_T profile to pT.dat\n" );

      file_pres = fopen( "pT_all.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         avPT[bin] /= nFrames;
         avPTc[bin] /= nFrames;
         avPT_K[bin] /= nFrames;
         avPT_U_LJ[bin] /= nFrames;
         avPT_U_LJc[bin] /= nFrames;
         avPT_U_QQ[bin] /= nFrames;
         avPT_U_14[bin] /= nFrames;
         avPT_U_BD1[bin] /= nFrames;
         avPT_U_BD2[bin] /= nFrames;
         fprintf ( file_pres, "%.3f  %.6f %.6f %.6f  %.6f %.6f %.6f  %.6f %.6f %.6f\n", 
                   avR[bin], avPT[bin], avPTc[bin], avPT_K[bin], 
                   avPT_U_LJ[bin], avPT_U_LJc[bin], avPT_U_QQ[bin],
                   avPT_U_14[bin], avPT_U_BD1[bin], avPT_U_BD2[bin] );
      }
      fclose( file_pres );

      file_pres = fopen( "pT_cut.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_pres, "%.3f  %.6f %.6f %.6f\n",
                   avR[bin], avPTc[bin], avPT_K[bin], avPTc[bin]-avPT_K[bin] );
      }
      fclose( file_pres );

      file_pres = fopen( "pT_full.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_pres, "%.3f  %.6f %.6f %.6f\n",
                   avR[bin], avPT[bin], avPT_K[bin], avPT[bin]-avPT_K[bin] );
      }
      fclose( file_pres );
/*
 *    write density.dat
 */
      printf( "   Printing density profile to density.dat\n" );

      file_dens = fopen( "density.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         avDens[bin] /= nFrames;
         fprintf ( file_dens, "%.3f  %.6f", avR[bin], avDens[bin] );
         for ( mol=0; mol<molTypes; mol++ )
         {
            avResDens[bin][mol] /= nFrames;
            fprintf ( file_dens, "  %.6f", avResDens[bin][mol] );
         }
         fprintf( file_dens, "\n" );
      }
      fclose( file_dens );


/*
 *    write dP.xvg
 */
      printf( "   Printing pN, pT, pN+r/2*dpN/dr to dP.xvg\n" );

      file_xvg = fopen( "dP.xvg", "w" );
      integr = 0.0;
      for ( bin=0; bin<maxBin; bin++ )
      {
/*
 *      Numerical differentiation
 *      http://en.wikipedia.org/wiki/Finite_difference_coefficients
 *                                 0     1     2     3     4     5     6
 *      Forward coefficients:   -49/20   6  -15/2  20/3 -15/4   6/5  -1/6
 *      Backward coefficients:   49/20  -6   15/2 -20/3  15/4  -6/5   1/6
 *                                -3    -2    -1     0     1     2     3
 *      Central coefficients:    -1/60  3/20 -3/4    0    3/4  -3/20  1/60
 */
         if ( bin < 3 )
         {
            avDPN[bin] = -49.0/20.0*avPN[bin]
                         + 6.0*avPN[bin+1] - 15.0/2.0*avPN[bin+2] + 20.0/3.0*avPN[bin+3]
                         - 15.0/4.0*avPN[bin+4] + 6.0/5.0*avPN[bin+5] - 1.0/6.0*avPN[bin+6];
            avDPN[bin] /= dR;
         }
         else if ( bin > maxBin-4 )
         {
            avDPN[bin] = 49.0/20.0*avPN[bin]
                         - 6.0*avPN[bin-1] + 15.0/2.0*avPN[bin-2] - 20.0/3.0*avPN[bin-3]
                         + 15.0/4.0*avPN[bin-4] - 6.0/5.0*avPN[bin-5] + 1.0/6.0*avPN[bin-6];
            avDPN[bin] /= dR;
         }
         else
         {
            avDPN[bin] = -1.0/60.0*avPN[bin-3] + 3.0/20.0*avPN[bin-2] - 3.0/4.0*avPN[bin-1]
                         + 3.0/4.0*avPN[bin+1] - 3.0/20.0*avPN[bin+2] + 1.0/60.0*avPN[bin+3];
            avDPN[bin] /= dR;
         }
         integr += avPT[bin] * avR[bin] * dR;
         fprintf ( file_xvg, "%.6f  %.6f  %.6f  %.6f  %.6f\n",
                   avR[bin], avPN[bin], avPT[bin], 
                   integr * 2.0 / pow(avR[bin],2.0),
                   avPN[bin]+0.5*avR[bin]*avDPN[bin] );
      }
      fclose( file_xvg );


/*
 *    fit density and compute surface tension
 */
      printf( "<> Fitting density profile and computing eff. surface tension...\n" );

      elg_param = elg_fit ( dR, maxBin, avR, avDens,
                            1.00, 0.01, boxSize*0.5-3.0, 0.2);

      dl = elg_param.liquid;
      dg = elg_param.gas;
      r0 = elg_param.r0;
      ksi = elg_param.ksi;

      dens_st = dens_to_st ( dR, maxBin, dl, dg, r0, ksi,
                             avR, avDens, avPN, avPNc );

      file_xvg = fopen( "fit_dens.xvg", "w" );
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_xvg, "%8.3f%20.10f%20.10f\n",
                   avR[bin], avDens[bin], 
                   dg+(dl-dg)*elg((avR[bin]-r0)/ksi) );
//                   0.5*(dl+dg)-0.5*(dl-dg)*elg((avR[bin]-r0)/ksi) );
      }
      fclose( file_xvg );

      file_xvg = fopen( "fit_densr2.xvg", "w" );
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_xvg, "%8.3f%20.10f%20.10f\n",
                   avR[bin], avDens[bin]*avR[bin]*avR[bin], 
                   (dg+(dl-dg)*elg((avR[bin]-r0)/ksi))*avR[bin]*avR[bin] );
//                   (0.5*(dl+dg)-0.5*(dl-dg)*elg((avR[bin]-r0)/ksi))*avR[bin]*avR[bin] );
      }
      fclose( file_xvg );

/*
 *    print results
 */
      printf( "   \n" );

      printf( "             LJ treatment  : %10s %10s            \n", "truncated", "full");
      printf( "   \n" );
      printf( "   Density of phase alpha  = %10.3f %10.3f  ( g/cm^3 )   \n", dens_st.dl, dens_st.dl );
      printf( "   Density of phase beta   = %10.3f %10.3f  ( g/cm^3 )   \n", dens_st.dg, dens_st.dg );
      printf( "   R of equimolar surface  = %10.2f %10.2f  ( nm     )   \n", dens_st.re, dens_st.re );
      printf( "   \n" );                                           

      printf( "   Work of formation       = %10.2f %10.2f  ( kJ/mol )  \n", 
                  dens_st.work_cut, 
                  dens_st.work_full );
      printf( "   Eff. surface tension    = %10.2f %10.2f  ( mJ/m^2 ) \n", 
                  dens_st.eff_st_cut / (nA*0.1), 
                  dens_st.eff_st_full / (nA*0.1) );
      printf( "   \n" );


/*
 *    fit pressure_N and compute surface tension
 */
      printf( "<> Fitting pressure_N profile and computing surface tension...\n" );

/*    
 *    full LJ interaction: 
 */
      elg_param = elg_fit ( dR, maxBin, avR, avPN,
                            30.0, 0.01, boxSize*0.5-3.0, 0.2);

      pl = elg_param.liquid;
      pg = elg_param.gas;
      r0 = elg_param.r0;
      ksi = elg_param.ksi;

      pres_st_full = pres_to_st ( dR, maxBin, pl, pg, r0, ksi,
                                  avR, avPN );

      file_xvg = fopen( "fit_pN.full.xvg", "w" );
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_xvg, "%8.3f%20.10f%20.10f\n",
                   avR[bin], avPN[bin],
                   pg+(pl-pg)*elg((avR[bin]-r0)/ksi) );
      }
      fclose( file_xvg );

      file_xvg = fopen( "fit_pNr2.full.xvg", "w" );
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_xvg, "%8.3f%20.10f%20.10f\n",
                   avR[bin], avPN[bin]*avR[bin]*avR[bin],
                   (pg+(pl-pg)*elg((avR[bin]-r0)/ksi))*avR[bin]*avR[bin] );
      }
      fclose( file_xvg );

/*    
 *    truncated LJ interaction: 
 */
      elg_param = elg_fit ( dR, maxBin, avR, avPNc,
                            30.0, 0.01, boxSize*0.5-3.0, 0.2);

      pl = elg_param.liquid;
      pg = elg_param.gas;
      r0 = elg_param.r0;
      ksi = elg_param.ksi;

      pres_st_cut = pres_to_st ( dR, maxBin, pl, pg, r0, ksi,
                                 avR, avPNc );

      file_xvg = fopen( "fit_pN.cut.xvg", "w" );
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_xvg, "%8.3f%20.10f%20.10f\n",
                   avR[bin], avPNc[bin],
                   pg+(pl-pg)*elg((avR[bin]-r0)/ksi) );
      }
      fclose( file_xvg );

      file_xvg = fopen( "fit_pNr2.cut.xvg", "w" );
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_xvg, "%8.3f%20.10f%20.10f\n",
                   avR[bin], avPNc[bin]*avR[bin]*avR[bin],
                   (pg+(pl-pg)*elg((avR[bin]-r0)/ksi))*avR[bin]*avR[bin] );
      }
      fclose( file_xvg );

/*
 *    print results
 */
      printf( "   \n" );

      printf( "             LJ treatment  : %10s %10s            \n", "truncated", "full");
      printf( "   \n" );

      printf( "   Pressure difference dP  = %10.2f %10.2f  ( bar    ) \n", 
                  ( pres_st_cut.pl - pres_st_cut.pg ) / (nA*0.01), 
                  ( pres_st_full.pl - pres_st_full.pg ) / (nA*0.01) );
      printf( "   Radius of Surface   Rs  = %10.2f %10.2f  ( nm     ) \n", 
                  pres_st_cut.rs, 
                  pres_st_full.rs );
      printf( "   Surface tension  gamma  = %10.2f %10.2f  ( mJ/m^2 ) \n", 
                  pres_st_cut.st / (nA*0.1),
                  pres_st_full.st / (nA*0.1) );
      printf( "   \n" );

/*
 *    end of main program
 */
      return 0;

   }

