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
         start_t = time(NULL);
         printf ( "\nPhase I: MPI execution\n" );
         printf ( "\n  Number of processors: %d\n", numprocs );
         printf ( "\n  Job started at %s", ctime(&start_t) );
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
         end_t = time(NULL);
         printf ( "\n  Job ended at %s", ctime(&end_t) );
         printf ( "\n  The MPI program terminated normally.\n");    

         delta_time = (int)(difftime(end_t,start_t));
         delta_hour = delta_time / 3600;
         delta_minute = (delta_time - delta_hour*3600) / 60;
         delta_second = delta_time - delta_hour*3600 - delta_minute*60;
         printf ( "\n  The calculation used %d hours %d minutes %d seconds.\n", 
                  delta_hour, delta_minute, delta_second );
      }

      if ( iproc == master )
      {
         printf ( "\nPhase II: Processing data\n" );
/*
 *       Run ./serial_process_data.x on master processor
 */ 
         sprintf ( command, "./serial_process_data.x %d > surftens.dat", numprocs );
         system ( command );
         printf ( "\nDone!\n" );
      }

      return 0;

   }

