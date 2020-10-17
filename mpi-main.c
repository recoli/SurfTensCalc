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
 *    MPI C program to call "sphere_calc_pres_dens.x"               *
 *    to compute surface tension of droplet                         *
 *    by Xin Li, TheoChem & Biology, KTH, Stockholm, Sweden         *
 *                   Version 2012-03-16                            *
 ********************************************************************/

   #include <stdio.h>
   #include <stdlib.h>
   #include <time.h>
   #include "mpi.h"

/********************************************************************
 *    Main program                                                  *
 ********************************************************************/

   int main( int argc, char *argv[] )
   {
/*
 *    There should be 5 arguments: nFrames, nStart, temperature, rCutoff, boxSize
 */
      if ( argc != 6 )
      {
         printf("%s\n", "Error: Incorrect number of arguments (should be 5)!");
         printf("%s\n", "1st argument: nFrame");
         printf("%s\n", "2nd argument: nStart");
         printf("%s\n", "3rd argument: temperature (in Kelvin)");
         printf("%s\n", "4th argument: cutoff radius (in nm)");
         printf("%s\n", "5th argument: approx. box size (in nm)");
         exit(1);
      }
      
      int nFrames, nStart;
      nFrames = atoi(argv[1]);
      nStart = atoi(argv[2]);
      
      double temperature;
      temperature = atof(argv[3]); /* K */

      double rCut;
      rCut = atof(argv[4]); /* nm */

      double boxSize;
      boxSize = atof(argv[5]); /* nm */

/*
 *    Define variables
 */
      char filename[20], command[256];
      int master = 0;
      int iproc, numprocs;
      MPI_Status status;
      time_t start_t, end_t;
      int delta_time, delta_hour, delta_minute, delta_second;
/*
 *    Initialize MPI
 */
      MPI_Init ( &argc, &argv );
/*
 *    Get the number of processes
 */
      MPI_Comm_size ( MPI_COMM_WORLD, &numprocs );
/*
 *    Determine the rank of this process
 */
      MPI_Comm_rank ( MPI_COMM_WORLD, &iproc );
/*
 *    Write to log
 */
      if ( iproc == master )
      {
         printf ( "*******************************************************************\n" );
         printf ( "*    C program to compute surface tension of spherical droplet    *\n" );
         printf ( "*                          by  Xin Li                             *\n" );
         printf ( "*           TheoChem & Biology, KTH, Stockholm, Sweden            *\n" );
         printf ( "*******************************************************************\n" );
         printf ( "\n" );
         printf ( "                       Version 2012-03-16                          \n" );
         printf ( "\n" );
         printf ( "Note:\n" );
         printf ( "a) Please check that the droplet does not cross box boundary during\n" );
         printf ( "   MD simulation.                                                  \n" );
         printf ( "b) It is recommended to use a long cutoff radius for Lennard-Jones \n" );
         printf ( "   interaction (e.g. 2.0 nm) during MD simulation, so that the long-\n" );
         printf ( "   range contribution to surface tension could be minimized.       \n" );
         printf ( "\n" );
         printf ( "Please cite the following papers:                                  \n" );
         printf ( "1) Thompson, S. M.; Gubbins, K. E.; Walton, J. P. R. B.; Chantry,  \n" );
         printf ( "   R. A. R. and Rowlinson, J. S.: A molecular dynamics study of    \n" );
         printf ( "   liquid drops, J. Chem. Phys., 81, 530-542, 1984.                \n" );
         printf ( "2) Li, X.; Hede, T.; Tu, Y.; Leck, C. and Agren, H.: Surface-active\n" );
         printf ( "   cis-pinonic acid in atmospheric droplets: A molecular dynamics  \n" );
         printf ( "   study, J. Phys. Chem. Lett., 1, 769-773, 2010.                  \n" );
         printf ( "3) Corti, D. S.; Kerr, K. J. and Torabi, K.: On the interfacial     \n" );
         printf ( "   thermodynamics of nanoscale droplets and bubbles, J. Chem. Phys.,\n" );
         printf ( "   135, 024701, 2011.                                               \n" );
         printf ( "4) Nakamura, T.; Shinoda, W. and Ikeshoji, T.: Novel numerical method\n" );
         printf ( "   for calculating the pressure tensor in spherical coordinates for  \n" );
         printf ( "   molecular systems, J. Chem. Phys., 135, 094106, 2011.             \n" );
         printf ( "\n" );

         start_t = time(NULL);
         printf ( "<> Job started at %s", ctime(&start_t) );
         printf ( "   %d frames will be read from the trajectory file.\n", nFrames );
         printf ( "   Only the last %d frames will be used in calculation.\n", nFrames-nStart );
         printf ( "   The temperature is %8.2f K\n", temperature );

         printf ( "<> Checking input files...\n" );
         system ( "./check_input.x > check_input.log" );
         system ( "cat check_input.log" );

         printf ( "<> Starting parallel calculations via MPI...\n" );
         printf ( "   Number of processors: %d\n", numprocs );
         printf ( "   Running sphere_calc_pres_dens.x on each processor...\n" );
      }
/*
 *    convert rank to filename
 */
      sprintf ( filename, "data-%d.log", iproc );
/*
 *    Run ./sphere_calc_pres_dens.x on each processor
 */
      sprintf ( command, "./sphere_calc_pres_dens.x %d %d %f %f %f > %s", 
                         nStart+(nFrames-nStart)/numprocs*(iproc+1), 
                         nStart+(nFrames-nStart)/numprocs*iproc, 
                         temperature, 
                         rCut,
                         boxSize,
                         filename );
      system ( command );
/*
 *    gzip log file
 */
      sprintf ( command, "gzip -f %s", filename );
      system ( command );
/*
 *    Shut down MPI
 */
      MPI_Finalize ( );
/*
 *    Write to log
 */
      if ( iproc == master )
      {
         end_t = time(NULL);
         printf ( "   Job ended at %s", ctime(&end_t) );

         delta_time = (int)(difftime(end_t,start_t));
         delta_hour = delta_time / 3600;
         delta_minute = (delta_time - delta_hour*3600) / 60;
         delta_second = delta_time - delta_hour*3600 - delta_minute*60;
         printf ( "   The calculation used %d hours %d minutes %d seconds.\n", 
                  delta_hour, delta_minute, delta_second );
      }

/*
 *    run ./sphere_process_data.x manually
 */ 

      return 0;

   }

