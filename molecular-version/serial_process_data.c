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
 *    There should be 2 argument: numprocs, temperature
 */
      if ( argc != 3 )
      {
         printf("%s\n", "Incorrect number of arguments (should be 2)!");
         printf("%s\n", "1st argument: numprocs");
         printf("%s\n", "2nd argument: temperature");
         exit(1);
      }
/*
 *    number of processors (data files)
 */
      int iproc, numprocs;
      numprocs = atoi(argv[1]);
/*
 *    temperature
 */
      double temperature;
      temperature = atof(argv[2]);
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
      const int maxBin = 1000;
      const double dR = 0.03; /*nm*/
      char filename[20], line[256];

      char text_work[256], text_pres[256], text_dens[256];
      FILE *file_data, *file_work, *file_pres, *file_dens;
      int read_text, read_data;

      double r, dens, pU, pN;
      double *aver_r, *aver_dens, *aver_pU, *aver_pN;

      int mol, molTypes;
      double *aver_molR, **aver_molDens;
      char *tmp;

      char command[256];
      FILE *file_fit_d, *file_fit_p;
      double fit_dl, fit_dg, fit_re, calc_work, eff_surf_tens;
/*
 *    set initial values and allocate arrays
 */
      molTypes = 20;

      aver_r = malloc( sizeof(double) * maxBin );
      aver_dens = malloc( sizeof(double) * maxBin );
      aver_pU = malloc( sizeof(double) * maxBin );
      aver_pN = malloc( sizeof(double) * maxBin );

      aver_molR = malloc(maxBin * sizeof(double));
      aver_molDens = malloc(maxBin * sizeof(double *));
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         aver_molDens[bin] = malloc(molTypes * sizeof(double));
      }

      for ( bin=0; bin<maxBin; bin++ )
      {
         aver_r[bin] = 0.0;
         aver_dens[bin] = 0.0;
         aver_pU[bin] = 0.0;
         aver_pN[bin] = 0.0;

         aver_molR[bin] = 0.0;
         for ( mol=0; mol<molTypes; mol++ )
         {
            aver_molDens[bin][mol] = 0.0;
         }
      }
/*
 *    read data files
 */
      sprintf ( text_work, "#Frame Work_of_formation\n" );
      sprintf ( text_pres, "#Radius Pressure_K Pressure_U Pressure_N\n" );
      sprintf ( text_dens, "#Radius Total_Density Individual_Densities\n" );

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
                  sscanf( line, "%lf%lf%lf%lf", &r, &dens, &pU, &pN );
                  aver_r[bin] += r;
                  aver_dens[bin] += dens;
                  aver_pU[bin] += pU;
                  aver_pN[bin] += pN;
               }
            }
            else if ( strcmp(line,text_dens) == 0 ) 
            {
               read_text = 3;
               read_data = 0;
               for ( bin=0; bin<maxBin; bin++ )
               {
                  fgets( line, sizeof(line), file_data );

                  tmp = strtok( line, " " );
                  aver_molR[bin] += atof(tmp);
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
      file_pres = fopen( "pressure.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         aver_r[bin] /= numprocs;
         aver_dens[bin] /= numprocs;
         aver_pU[bin] /= numprocs;
         aver_pN[bin] /= numprocs;
         fprintf ( file_pres, "%8.3f%20.10f%20.10f%20.10f\n", 
                   aver_r[bin], aver_dens[bin], aver_pU[bin], aver_pN[bin] );
      }
      fclose( file_pres );
/*
 *    write density.dat
 */
      file_dens = fopen( "density.dat", "w" ) ;
      for ( bin=0; bin<maxBin; bin++ )
      {
         aver_molR[bin] /= numprocs;
         fprintf ( file_dens, "%8.3f", aver_molR[bin] );
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
      sprintf ( command, "./fit_dens.x > fit_dens.log" );
      system ( command );
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

      calc_work = 0.0;
      for ( bin=0; bin<maxBin; bin++ )
      {
         calc_work += ( aver_pN[bin]*aver_r[bin]*aver_r[bin] ) * dR;
      }
      calc_work *= pi * 2.0 ;  /*kJ/mol*/
      eff_surf_tens = calc_work * 3.0 / ( 4.0 * pi * fit_re * fit_re ) ; /*kJ/mol/nm^2*/
      eff_surf_tens /= (nA*0.1) ;  /*mN/m*/
      calc_work /= (nA*10.0) ;  /*10^-19 J*/

      printf( "   Density of the cluster     = %f nm^-3\n",   fit_dl );
      printf( "   Equimolar dividing surface = %f nm\n",      fit_re );
      printf( "   Work of formation          = %f e-19 J\n",  calc_work );
      printf( "   Effective surface tension  = %f mJ m^-2\n", eff_surf_tens );

/*
 *    end of main program
 */
      return 0;

   }

