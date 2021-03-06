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
 *    MPI C program to call "serial_calc_pres_dens.x"               *
 *    to compute surface tension of droplet                         *
 *    by Xin Li, TheoChem & Biology, KTH, Stockholm, Sweden         *
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
 *    There should be 3 arguments: nFrames, nStart, temperature
 */
      if ( argc != 4 )
      {
         printf("%s\n", "Incorrect number of arguments (should be 3)!");
         printf("%s\n", "1st argument: nFrame");
         printf("%s\n", "2nd argument: nStart");
         printf("%s\n", "3rd argument: temperature (in Kelvin)");
         exit(1);
      }
      
      int nFrames, nStart;
      nFrames = atoi(argv[1]);
      nStart = atoi(argv[2]);
      
      double temperature;
      temperature = atof(argv[3]); /* K */
/*
 *    Define variables
 */
      char filename[10], command[256];
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
         printf ( "*****************************************************\n" );
         printf ( "*  C program to compute surface tension of droplet  *\n" );
         printf ( "*           Xin Li, TheoChem & Biology, KTH         *\n" );
         printf ( "*               KTH, Stockholm, Sweden              *\n" );
         printf ( "*****************************************************\n" );

         printf ( "\nThis program is based on the following papers:\n" );
         printf ( "   Thompson, et al., J. Chem. Phys., 1984, 81, 530-542.\n" );
         printf ( "   Brodskaya, et al., J. Colloid Interface Sci., 1996, 180, 86-97.\n" );
         printf ( "   Li, et al., J. Phys. Chem. Lett. 2010, 1, 769-773.\n" );

         start_t = time(NULL);
         printf ( "\nStep I: MPI execution\n" );
         printf ( "   There are %d frames; only the last %d frames will be used.\n", nFrames, nFrames-nStart );
         printf ( "   The temperature is %8.2f K\n", temperature );
         printf ( "<> Checking input files...\n" );
         fflush(stdout);

         system ( "./check_input.x > check_input.log" );
         system ( "cat check_input.log" );

         printf ( "<> Job started at %s", ctime(&start_t) );
         printf ( "   Number of processors: %d\n", numprocs );
         printf ( "   Running serial_calc_pres_dens.x on each processor...\n" );
         fflush(stdout);
      }
/*
 *    convert rank to filename
 */
      sprintf ( filename, "%d", iproc );
/*
 *    Run ./serial_calc_pres_dens.x on each processor
 */
      sprintf ( command, "./serial_calc_pres_dens.x %d %d %f > data-%s.log", 
                         nStart+(nFrames-nStart)/numprocs*(iproc+1), 
                         nStart+(nFrames-nStart)/numprocs*iproc, 
                         temperature, 
                         filename );
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
         printf ( "\nStep II: Processing data\n" );
         fflush(stdout);
/*
 *       Run ./serial_process_data.x on master processor
 */ 
         sprintf ( command, "./serial_process_data.x %d %f > surftens.dat", 
                            numprocs, temperature );
         system ( command );
         system ( "cat surftens.dat" );

         end_t = time(NULL);
         printf ( "\n<> Job ended at %s", ctime(&end_t) );

         delta_time = (int)(difftime(end_t,start_t));
         delta_hour = delta_time / 3600;
         delta_minute = (delta_time - delta_hour*3600) / 60;
         delta_second = delta_time - delta_hour*3600 - delta_minute*60;
         printf ( "The calculation used %d hours %d minutes %d seconds.\n", 
                  delta_hour, delta_minute, delta_second );
         fflush(stdout);
      }

      return 0;

   }

