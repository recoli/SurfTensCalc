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
 *    Serial C program to fit density profile to tanh function      *
 *    by Xin Li, TheoChem & Biology, KTH, Stockholm, Sweden         *
 ********************************************************************/

   #include <math.h>
   #include <stdio.h>
   #include <stdlib.h>
   #include <string.h>
   #include <stdbool.h>

/********************************************************************
 *    Main program                                                  *
 ********************************************************************/

   int main( int argc, char *argv[] )
   {

/********************************************************************
 *    Define variables                                              *
 ********************************************************************/

      int bin;
      const int  maxbin = 1000;
      const double delr = 0.03;
      
      double *r, *dens, *dd;
      double dl, dg, integr, re;

      FILE   *file_dat, *file_xvg ;
      char   line[256];
 
      int    count;

/********* fit the density profile to hyperbolic tangent function ************
 *  
 *    dens(r) = 0.5 * ( dens_l + dens_g ) - 
 *              0.5 * ( dens_l - dens_g ) * tanh( ( r - r0 ) / ksi )
 *  
 *    beta(1) = dens_l       :  density of the liquid phase
 *    beta(2) = dens_g       :  density of the gas phase
 *    beta(3) = r0           :  radius of the drop
 *    beta(4) = ksi          :  thickness parameter of the interface
 *    iter                   :  iterate or not
 *  
 *****************************************************************************/

      bool         iter;
      double       beta[4];
      double       guess_r;    /* initial guess of droplet radius*/
      const double damp = 0.8; /* damping factor in iteration*/
      double       *y;
      double       **jaco;

      int          i, j, k;
      double       a[4][4], b[4];
      double       aii, aji;

/*
 *    allocate storage for arrays
 */
      r    = malloc(maxbin * sizeof(double));
      dens = malloc(maxbin * sizeof(double));
      dd   = malloc(maxbin * sizeof(double));
      y    = malloc(maxbin * sizeof(double));
/*
 *    allocate storage for an array of pointers
 */
      jaco = malloc(maxbin * sizeof(double *));
/*
 *    for each pointer, allocate storage for an array of ints
 */
      for (bin=0; bin<maxbin; bin++) {
         jaco[bin] = malloc(4 * sizeof(double));
      }

/*
 *    read density profile   
 */ 
      for ( bin=0; bin<maxbin; bin++ ) 
      {
         r[bin] = 0.0;
         dens[bin] = 0.0;
      }

      file_dat = fopen("density.dat","r") ;
      if ( NULL==file_dat ) 
      {
         printf( "Cannot open file: density.dat !\n" ) ;
         exit(1);
      }
      for ( bin=0; bin<maxbin; bin++ ) 
      {
         fgets( line, sizeof( line ), file_dat );
         sscanf( line, "%lf%lf", &r[bin], &dens[bin] );
      }
      fclose( file_dat );
/*
 *    guess the radius of droplet
 */
      for ( bin=0; bin<maxbin; bin++ ) 
      {
         if ( dens[bin]>5.0 && dens[bin]<10.0 && bin>20 )
         {
            guess_r = r[bin];
            break;
         }
      }
/*
 *    initialize beta vector  
 */
      beta[0] =   33.50;
      beta[1] =    0.01;
      beta[2] = guess_r;
      beta[3] =    0.30; 

/*
 *    iteratively solve the non-linear least square fitting  
 */ 
      count = 0;
      iter = true;
      while (iter) {
/*
 *       calculate y(bin) and jacobian matrices 
 *       dens(r) = 0.5 * ( dens_l + dens_g ) 
 *               - 0.5 * ( dens_l - dens_g ) * tanh( ( r - r0 ) / ksi )  
 *
 *       beta[0] ---> dens_l
 *       beta[1] ---> dens_g
 *       beta[2] ---> r0
 *       beta[3] ---> ksi
 */
         for ( bin=0; bin<maxbin; bin++ ) {
            y[bin]       = 0.5 * ( beta[0] + beta[1] )  
                         - 0.5 * ( beta[0] - beta[1] )  
                         * tanh( ( r[bin]-beta[2] ) / beta[3] );
            jaco[bin][0] = 0.5 - 0.5 * tanh( ( r[bin]-beta[2] ) / beta[3] );
            jaco[bin][1] = 0.5 + 0.5 * tanh( ( r[bin]-beta[2] ) / beta[3] );
            jaco[bin][2] = -0.5 * ( beta[0] - beta[1] )  
                         / pow( cosh( ( r[bin]-beta[2] )/beta[3] ), 2.0 )
                         * ( -1.0 / beta[3] );
            jaco[bin][3] = -0.5 * ( beta[0] - beta[1] )  
                         / pow( cosh( ( r[bin]-beta[2] )/beta[3] ), 2.0 ) 
                         * ( r[bin]-beta[2] ) * ( -1.0/pow(beta[3],2.0) );
         }
/*
 *       determine coefficients of the linear equation system  
 */
         for ( i=0; i<4; i++ ) 
         {
            for ( j=0; j<4; j++ ) 
            {
               a[i][j] = 0.0;
            }
            b[i] = 0.0;
         }
/*
 *       summation starts from bin=5 to avoid bad statistics near the center 
 */
         for ( bin=5; bin<maxbin; bin++ ) 
         {
            for ( i=0; i<4; i++ ) 
            {
               for ( j=0; j<4; j++ ) 
               {
                  a[i][j] += jaco[bin][i] * jaco[bin][j];
               }
               b[i] += jaco[bin][i] * ( dens[bin] - y[bin] );
            }
         }
/*
 *       use Gaussian elimination  
 */
         for ( i=0; i<4; i++ ) 
         {
            aii = a[i][i];
            for ( j=0; j<4; j++ ) 
            {
               a[i][j] /= aii;
            }
            b[i] /= aii;
            for ( j=i+1; j<4; j++ ) 
            {
               aji = a[j][i];
               for ( k=0; k<4; k++ ) 
               {
                  a[j][k] -= aji * a[i][k];
               }
               b[j] -= aji * b[i];
            }
         }
         for ( i=4-1; i>=0; i-- ) 
         {
            b[i] /= a[i][i];
            for ( j=i-1; j>=0; j-- ) 
            {
               b[j] -= a[j][i] * b[i];
            }
         }
/*
 *       generate new beta vector and check the change  
 */ 
         for ( i=0; i<4; i++ ) 
         {
            beta[i] += b[i] * damp;
            if ( fabs( b[i]/beta[i] ) > 1.0e-8 ) 
            {
               iter = true;
               break;
            } 
            else 
            {
               iter = false;
               continue;
            }
         }

         count++;
         if ( count > 5000 ) 
         {
            printf("Error in fit_dens: Iteration exceeds 5000!\n");
            exit(1);
         }
/*
 *    end of loop (iter)  
 */
      }

/*
 *    compare the calculated and the fitted density profile  
 */
      file_xvg = fopen( "fit_dens.xvg", "w" );
      for ( bin=0; bin<maxbin; bin++ ) 
      {
         y[bin] = 0.5*(beta[0]+beta[1]) - 0.5*(beta[0]-beta[1]) 
                * tanh((r[bin]-beta[2])/beta[3]);
         fprintf ( file_xvg, "%8.3f%20.10f%20.10f\n", r[bin], dens[bin], y[bin] );
      }
      fclose( file_xvg );

/*
 *    calculate radius of equimolar dividing surface  
 */
      dl = beta[0];
      dg = beta[1];
      integr = 0.0;
      for ( bin=0; bin<maxbin; bin++ ) 
      {
/*
 *       numerical differentiation
 *                -f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)
 *       f'(x) = ----------------------------------------
 *                                 12h                       
 */
         if ( bin==0 ) 
         {
            dd[bin] = ( dens[bin+1] - dens[bin] ) / delr;
         } 
         else if ( bin==maxbin-1 ) 
         {
            dd[bin] = ( dens[bin] - dens[bin-1] ) / delr;
         } 
         else if ( bin==1 || bin==maxbin-2 ) 
         {
            dd[bin] = ( dens[bin+1] - dens[bin-1] ) / ( 2.0 * delr );
         } 
         else 
         {
            dd[bin] = -dens[bin+2] + 8.0 * dens[bin+1] 
                    - 8.0 * dens[bin-1] + dens[bin-2];
            dd[bin] /= ( 12.0 * delr );
         }
         integr += pow(r[bin],3.0) * dd[bin] * delr;
      }
      re = integr / ( dg - dl );
      re = pow(re,0.3333333333); /* nm */

/*
 *    write final results
 */
      printf( "%20.10f  (%-s)\n", dl, "Liquid_Density, nm^-3" ) ;
      printf( "%20.10f  (%-s)\n", dg, "Gas_Density, nm^-3" ) ;
      printf( "%20.10f  (%-s)\n", re, "Re, nm" ) ;

/*
 *    release arrays
 */
      free(r);
      free(dens);
      free(dd);
      free(y);
/*
 *    release each pointer in an array
 */
      for (bin=0; bin<maxbin; bin++) 
      {
         free(jaco[bin]);
      }
/*
 *    release 2D array
 */
      free(jaco);

/*****************************************************************************
 *    End of Program                                                         *
 *****************************************************************************/

      return 0;
      }

