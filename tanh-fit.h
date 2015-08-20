/*****************************************************************************
 This file is part of the SurfTensCalc program.                                      
 Copyright (C) 2012 Xin Li
                                                                           
 Filename:  tanh-fit.h
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
 *****************************************************************************/


/********* fit the pressure_N profile to hyperbolic tangent function *********
 *                                                                           *
 *    pres(r) = 0.5 * ( pres_l + pres_g ) -                                  *
 *              0.5 * ( pres_l - pres_g ) * tanh( ( r - r0 ) / ksi )         *
 *                                                                           *
 *    beta(1) = pres_l       :  pressure of the liquid phase                 *
 *    beta(2) = pres_g       :  pressure of the gas phase                    *
 *    beta(3) = r0           :  radius of the drop                           *
 *    beta(4) = ksi          :  thickness parameter of the interface         *
 *    iter                   :  iterate or not                               *
 *                                                                           *
 *****************************************************************************/

   typedef struct {
      double liquid, gas, r0, ksi;
   } Vec_tanh;

   Vec_tanh tanh_fit ( double delr, int maxbin, double r[], double pn[],
                       double b_0, double b_1, double b_2, double b_3 )
   {
      Vec_tanh tanh_param;
/*
 *    damping factor for iteration
 */
      const double damp = 0.8; 
/*
 *    variables
 */
      int bin, count;
      int iter;
      double beta[4], *y, **jaco;

      int i, j, k;
      double a[4][4], b[4], aii, aji;
/*
 *    allocate arrays
 */
      y   = malloc(maxbin * sizeof(double));
      jaco = malloc(maxbin * sizeof(double *));
      for (bin=0; bin<maxbin; bin++) {
         jaco[bin] = malloc(4 * sizeof(double));
      }
/*
 *    initial guess for beta vector  
 */
      beta[0] = b_0;
      beta[1] = b_1;
      beta[2] = b_2;
      beta[3] = b_3;
/*
 *    iteratively solve the non-linear least square fitting  
 */ 
      count = 0;
      iter = 1;
      while (iter==1) {
/*
 *       calculate y(bin) and jacobian matrices 
 *       pres(r) = 0.5 * ( pres_l + pres_g ) 
 *               - 0.5 * ( pres_l - pres_g ) * tanh( ( r - r0 ) / ksi )  
 *
 *       beta[0] ---> pres_l
 *       beta[1] ---> pres_g
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
               b[i] += jaco[bin][i] * ( pn[bin] - y[bin] );
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
               iter = 1;
               break;
            } 
            else 
            {
               iter = 0;
               continue;
            }
         }

         count++;
         if ( count > 5000 ) 
         {
            printf("<> Error: Iteration exceeds 5000 in tanh_fit !\n");
            exit(1);
         }
/*
 *    end of loop (iter)  
 */
      }
/*
 *    free arrays
 */
      free(y);
      free(jaco);

      tanh_param.liquid = beta[0];
      tanh_param.gas = beta[1];
      tanh_param.r0 = beta[2];
      tanh_param.ksi = beta[3];

      return tanh_param;
      }


/******* fit the pressure_N*r^2 profile to hyperbolic tangent function *******
 *                                                                           *
 *    pres(r) = ( 0.5 * ( pres_l + pres_g ) -                                *
 *                0.5 * ( pres_l - pres_g ) * tanh( ( r - r0 ) / ksi ) )     *
 *              * r^2                                                        *
 *                                                                           *
 *    beta(1) = pres_l       :  pressure of the liquid phase                 *
 *    beta(2) = pres_g       :  pressure of the gas phase                    *
 *    beta(3) = r0           :  radius of the drop                           *
 *    beta(4) = ksi          :  thickness parameter of the interface         *
 *    iter                   :  iterate or not                               *
 *                                                                           *
 *****************************************************************************/

   Vec_tanh tanh_r2_fit ( double delr, int maxbin, double r[], double pn[],
                          double b_0, double b_1, double b_2, double b_3 )
   {
      Vec_tanh tanh_param;
/*
 *    damping factor for iteration
 */
      const double damp = 0.8; 
/*
 *    variables
 */
      int bin, count;
      int iter;
      double beta[4], *y, **jaco;

      int i, j, k;
      double a[4][4], b[4], aii, aji;
/*
 *    allocate arrays
 */
      y   = malloc(maxbin * sizeof(double));
      jaco = malloc(maxbin * sizeof(double *));
      for (bin=0; bin<maxbin; bin++) {
         jaco[bin] = malloc(4 * sizeof(double));
      }
/*
 *    initial guess for beta vector  
 */
      beta[0] = b_0;
      beta[1] = b_1;
      beta[2] = b_2;
      beta[3] = b_3;
/*
 *    iteratively solve the non-linear least square fitting  
 */ 
      count = 0;
      iter = 1;
      while (iter==1) {
/*
 *       calculate y(bin) and jacobian matrices 
 *       pres(r) = ( 0.5 * ( pres_l + pres_g ) 
 *                 - 0.5 * ( pres_l - pres_g ) * tanh( ( r - r0 ) / ksi ) ) * r^2
 *
 *       beta[0] ---> pres_l
 *       beta[1] ---> pres_g
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

            y[bin] *= r[bin]*r[bin];
            jaco[bin][0] *= r[bin]*r[bin];
            jaco[bin][1] *= r[bin]*r[bin];
            jaco[bin][2] *= r[bin]*r[bin];
            jaco[bin][3] *= r[bin]*r[bin];
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
 *       summation starts from bin=0 
 *       no need to avoid bad statistics near the center 
 */
         for ( bin=0; bin<maxbin; bin++ ) 
         {
            for ( i=0; i<4; i++ ) 
            {
               for ( j=0; j<4; j++ ) 
               {
                  a[i][j] += jaco[bin][i] * jaco[bin][j];
               }
               b[i] += jaco[bin][i] * ( pn[bin]*r[bin]*r[bin] - y[bin] );
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
               iter = 1;
               break;
            } 
            else 
            {
               iter = 0;
               continue;
            }
         }

         count++;
         if ( count > 5000 ) 
         {
            printf("<> Error: Iteration exceeds 5000 in tanh_r2_fit !\n");
            exit(1);
         }
/*
 *    end of loop (iter)  
 */
      }
/*
 *    free arrays
 */
      free(y);
      free(jaco);

      tanh_param.liquid = beta[0];
      tanh_param.gas = beta[1];
      tanh_param.r0 = beta[2];
      tanh_param.ksi = beta[3];

      return tanh_param;
      }

