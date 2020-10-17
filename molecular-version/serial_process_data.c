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
 *    C program to analyze data from serial_calc_pres_dens.x        *
 *    by Xin Li, TheoChem & Biology, KTH, Stockholm, Sweden         *
 ********************************************************************/

   #include <math.h>
   #include <stdio.h>
   #include <stdlib.h>
   #include <string.h>

/********************************************************************
 *    Main program                                                  *
 ********************************************************************/

   int main( int argc, char *argv[] )
   {
/*
 *    There should be 1 argument: numprocs
 */
      if ( argc != 3 )
      {
         printf("%s\n", "Incorrect number of arguments (should be 2)!");
         printf("%s\n", "1st argument: numprocs");
         printf("%s\n", "2nd argument: approx. box size (in nm)");
         exit(1);
      }
/*
 *    number of processors (data files)
 */
      int iproc, numprocs;
      numprocs = atoi(argv[1]);

      double boxSize;
      boxSize = atof(argv[2]); /* nm */

/*
 *    Define variables
 */

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
      const double dR = 0.03; /*nm*/
      int maxBin ;
      maxBin = (int)( boxSize * 0.4 / dR ) ;

      char filename[20], line[256];

      char text_work[256], text_pres[256], text_dens[256];
      FILE *file_data, *file_work, *file_pres, *file_dens;
      int read_text, read_data;

/*
      double r, pK, pU, pN;
      double *aver_r, *aver_pK, *aver_pU, *aver_pN;
*/
      double r, pN_full, pN_cut, pK, pULJ, pUQQ, pULJc;
      double *aver_r, *aver_pN_full, *aver_pN_cut;
      double *aver_pK, *aver_pULJ, *aver_pUQQ, *aver_pULJc;

      int mol, dens, molTypes;
      double *aver_molR, *aver_dens, **aver_molDens;
      char *tmp;

      char command[256];
      FILE *file_fit_d, *file_fit_p;

      double fit_dl, fit_dg, fit_re;

      double calc_work_cut, eff_surf_tens_cut;
      double fit_pl_cut, fit_pg_cut, fit_rs_cut, fit_st_cut, fit_work_cut, deltaP_cut;

      double calc_work_full, eff_surf_tens_full;
      double fit_pl_full, fit_pg_full, fit_rs_full, fit_st_full, fit_work_full, deltaP_full;

/*
 *    set initial values and allocate arrays
 */
      molTypes = 20;

      aver_r = malloc( sizeof(double) * maxBin );
      aver_pN_full = malloc( sizeof(double) * maxBin );
      aver_pN_cut = malloc( sizeof(double) * maxBin );
      aver_pK = malloc( sizeof(double) * maxBin );
      aver_pULJ = malloc( sizeof(double) * maxBin );
      aver_pUQQ = malloc( sizeof(double) * maxBin );
      aver_pULJc = malloc( sizeof(double) * maxBin );

      aver_molR = malloc(maxBin * sizeof(double));
      aver_dens = malloc( sizeof(double) * maxBin );
      aver_molDens = malloc(maxBin * sizeof(double *));
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         aver_molDens[bin] = malloc(molTypes * sizeof(double));
      }

      for ( bin=0; bin<maxBin; bin++ )
      {
         aver_r[bin] = 0.0;
         aver_pN_full[bin] = 0.0;
         aver_pN_cut[bin] = 0.0;
         aver_pK[bin] = 0.0;
         aver_pULJ[bin] = 0.0;
         aver_pUQQ[bin] = 0.0;
         aver_pULJc[bin] = 0.0;

         aver_molR[bin] = 0.0;
         aver_dens[bin] = 0.0;
         for ( mol=0; mol<molTypes; mol++ )
         {
            aver_molDens[bin][mol] = 0.0;
         }
      }
/*
 *    read data files
 */

      printf( "<> Collecting data from each processor...\n" );

      sprintf ( text_work, "#Frame  W_Full  W_Cut\n" );
      sprintf ( text_pres, "#Radius  P_N_full  P_N_cut  P_K  P_U_LJ  P_U_QQ  P_U_LJc\n" );
      sprintf ( text_dens, "#Radius Total_Density Individual_Densities\n" );

      printf( "   Printing work of formation to work.dat\n" );

      file_work = fopen( "work.dat", "w" ) ;

      for ( iproc=0; iproc<numprocs; iproc++ )
      {
         sprintf ( filename, "data-%d.log", iproc );
         file_data = fopen( filename, "r" ) ;
         if ( NULL==file_data ) 
         {
            printf( "Cannot open file: %s !\n", filename ) ;
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
            }
            else if ( strcmp(line,text_pres) == 0 ) 
            {
               read_text = 2;
               read_data = 0;
               for ( bin=0; bin<maxBin; bin++ )
               {
                  fgets( line, sizeof(line), file_data );
                  sscanf( line, "%lf%lf%lf%lf%lf%lf%lf", 
                          &r, &pN_full, &pN_cut, &pK, &pULJ, &pUQQ, &pULJc );
                  aver_r[bin] += r;
                  aver_pN_full[bin] += pN_full;
                  aver_pN_cut[bin] += pN_cut;
                  aver_pK[bin] += pK;
                  aver_pULJ[bin] += pULJ;
                  aver_pUQQ[bin] += pUQQ;
                  aver_pULJc[bin] += pULJc;
               }
            }
            else if ( strcmp(line,text_dens) == 0 ) 
            {
               read_text = 3;
               read_data = 0;
               for ( bin=0; bin<maxBin; bin++ )
               {
                  fgets( line, sizeof(line), file_data );
/*
 *                read the 1st column: radius
 */
                  tmp = strtok( line, " " );
                  aver_molR[bin] += atof(tmp);
/*
 *                read the 2nd column: total density
 */
                  tmp = strtok( NULL, " " );
                  if ( NULL==tmp )
                  {
                     break;
                  }
                  aver_dens[bin] += atof(tmp);
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
                     aver_molDens[bin][mol] += atof(tmp);
                     mol ++;
                  } 
                  molTypes = mol;

               }
            }
            else
            {
               read_data = 1;
            }

            if ( read_text==1 && read_data==1 ) 
            {
               fprintf( file_work, "%s", line );
            }

         }

         fclose( file_data );

      }
/*
 *    write pressure.dat
 */

      printf( "   Printing pressure profile to pressure.dat\n" );

      file_pres = fopen( "pressure_all.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         aver_r[bin] /= numprocs;
         aver_pN_full[bin] /= numprocs;
         aver_pN_cut[bin] /= numprocs;
         aver_pK[bin] /= numprocs;
         aver_pULJ[bin] /= numprocs;
         aver_pUQQ[bin] /= numprocs;
         aver_pULJc[bin] /= numprocs;
         fprintf ( file_pres, "%8.3f%20.10f%20.10f%20.10f%20.10f%20.10f%20.10f\n", 
                   aver_r[bin], aver_pN_full[bin], aver_pN_cut[bin], 
                   aver_pK[bin], aver_pULJ[bin], aver_pUQQ[bin], aver_pULJc[bin] );
      }
      fclose( file_pres );

      file_pres = fopen( "pressure_cut.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_pres, "%8.3f%20.10f%20.10f%20.10f\n",
                   aver_r[bin], aver_pK[bin], aver_pUQQ[bin]+aver_pULJc[bin], aver_pN_cut[bin] );
      }
      fclose( file_pres );

      file_pres = fopen( "pressure_full.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         fprintf ( file_pres, "%8.3f%20.10f%20.10f%20.10f\n",
                   aver_r[bin], aver_pK[bin], aver_pUQQ[bin]+aver_pULJ[bin], aver_pN_full[bin] );
      }
      fclose( file_pres );

/*
 *    write density.dat
 */

      printf( "   Printing density profile to density.dat\n" );

      file_dens = fopen( "density.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         aver_molR[bin] /= numprocs;
         aver_dens[bin] /= numprocs;
         fprintf ( file_dens, "%8.3f%20.10f", aver_molR[bin], aver_dens[bin] );
         for ( mol=0; mol<molTypes; mol++ )
         {
            aver_molDens[bin][mol] /= numprocs;
            fprintf ( file_dens, "%20.10f", aver_molDens[bin][mol] );
         }
         fprintf( file_dens, "\n" );
      }
      fclose( file_dens );

/*
 *    fit density and compute surface tension
 */

      printf( "<> Fitting density profile and computing surface tension...\n" );

      system( "./fit_dens.x > fit_dens.log" );
      file_fit_d = fopen( "fit_dens.log", "r" ) ;
      if ( NULL==file_fit_d )
      {
         printf( "Cannot open file: fit_dens.log!\n" ) ;
         exit(1);
      }
      fgets( line, sizeof(line), file_fit_d );
      sscanf( line, "%lf", &fit_dl );
      fgets( line, sizeof(line), file_fit_d );
      sscanf( line, "%lf", &fit_dg );
      fgets( line, sizeof(line), file_fit_d );
      sscanf( line, "%lf", &fit_re );
      fclose( file_fit_d );


/*
 *    fit pressure_N and compute surface tension
 */

      printf( "<> Fitting pressure_N profile and computing surface tension...\n" );

/*    truncated LJ interaction: */

      system( "cp pressure_cut.dat pressure.dat" );

      sprintf ( command, "./fit_pN.x %f > fit_pN.log", boxSize );
      system ( command );

      file_fit_p = fopen( "fit_pN.log", "r" ) ;
      if ( NULL==file_fit_p )
      {
         printf( "Cannot open file: fit_pN.log!\n" ) ;
         exit(1);
      }
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_pl_cut );
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_pg_cut );
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_rs_cut );
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_st_cut );
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_work_cut );

      deltaP_cut = ( fit_pl_cut - fit_pg_cut ) / nA; /* e+07 Pa */

/*    full LJ interaction: */

      system( "cp pressure_full.dat pressure.dat" );

      sprintf ( command, "./fit_pN.x %f > fit_pN.log", boxSize );
      system ( command );

      file_fit_p = fopen( "fit_pN.log", "r" ) ;
      if ( NULL==file_fit_p )
      {
         printf( "Cannot open file: fit_pN.log!\n" ) ;
         exit(1);
      }
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_pl_full );
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_pg_full );
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_rs_full );
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_st_full );
      fgets( line, sizeof(line), file_fit_p );
      sscanf( line, "%lf", &fit_work_full );

      deltaP_full = ( fit_pl_full - fit_pg_full ) / nA; /* e+07 Pa */

/*
 *    compute work of formation and effective surface tension
 */

      calc_work_cut = 0.0;
      calc_work_full = 0.0;
      for ( bin=0; bin<maxBin; bin++ )
      {
         calc_work_cut += ( (aver_pN_cut[bin]-fit_pg_cut)*aver_r[bin]*aver_r[bin] ) * dR;
         calc_work_full += ( (aver_pN_full[bin]-fit_pg_full)*aver_r[bin]*aver_r[bin] ) * dR;
      }

      calc_work_cut *= pi * 2.0 ;  /*kJ/mol*/
      eff_surf_tens_cut = calc_work_cut * 3.0 / ( 4.0 * pi * fit_re * fit_re ) ; /*kJ/mol/nm^2*/
      eff_surf_tens_cut /= (nA*0.1) ;  /*mN/m*/
      calc_work_cut /= (nA*10.0) ;  /*10^-19 J*/

      calc_work_full *= pi * 2.0 ;  /*kJ/mol*/
      eff_surf_tens_full = calc_work_full * 3.0 / ( 4.0 * pi * fit_re * fit_re ) ; /*kJ/mol/nm^2*/
      eff_surf_tens_full /= (nA*0.1) ;  /*mN/m*/
      calc_work_full /= (nA*10.0) ;  /*10^-19 J*/

      printf( "   \n" );
      printf( "   (Hyperbolic tangent fit of density)\n" );
      printf( "   Density of the cluster     = %f nm^-3\n",   fit_dl );
      printf( "   Equimolar dividing surface = %f nm\n",      fit_re );
      printf( "   \n" );
      printf( "   (Hyperbolic tangent fit of pressure_N)\n" );
      printf( "   Pressure difference dP     = %f e+07 Pa (%f)\n", deltaP_cut, deltaP_full );
      printf( "   (Mechanical route)\n" );
      printf( "   Surface tension  gamma     = %f mJ m^-2 (%f)\n", fit_st_cut, fit_st_full );
      printf( "   (Young-Laplace equation)\n" );
      printf( "   Radius of Surface   Rs     = %f nm (%f)\n",      fit_rs_cut, fit_rs_full );
      printf( "   \n" );
      printf( "   (W = 4pi/3 * gamma * Rs^2)\n" );
      printf( "   Work of formation    W     = %f e-19 J (%f)\n",  fit_work_cut, fit_work_full );
      printf( "   \n" );
      printf( "   (For reference)\n" );
      printf( "   Work of formation (ref)    = %f e-19 J (%f)\n",  calc_work_cut, calc_work_full );
      printf( "   Eff. surface tension (ref) = %f mJ m^-2 (%f)\n", eff_surf_tens_cut, eff_surf_tens_full );
      printf( "   \n" );


/*
 *    end of main program
 */
      return 0;

   }

