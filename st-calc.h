/*****************************************************************************
 This file is part of the SurfTensCalc program.                                      
 Copyright (C) 2012 Xin Li
                                                                           
 Filename:  st-calc.h
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


   const double pi   =   3.14159265358979323846;
   const double nA   =   6.02214129 ;  /* 10^23 mol^-1 */


   double integr_r3 ( double dR, int maxBin, double r[], double pn[] )
   {
      double *dpn;
      dpn = malloc(maxBin * sizeof(double));
 
      double integr;
      integr = 0.0;

      int bin;

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
            dpn[bin] = -49.0/20.0*pn[bin]
                       + 6.0*pn[bin+1] - 15.0/2.0*pn[bin+2] + 20.0/3.0*pn[bin+3]
                       - 15.0/4.0*pn[bin+4] + 6.0/5.0*pn[bin+5] - 1.0/6.0*pn[bin+6];
            dpn[bin] /= dR;
         }
         else if ( bin > maxBin-4 )
         {
            dpn[bin] = 49.0/20.0*pn[bin]
                       - 6.0*pn[bin-1] + 15.0/2.0*pn[bin-2] - 20.0/3.0*pn[bin-3]
                       + 15.0/4.0*pn[bin-4] - 6.0/5.0*pn[bin-5] + 1.0/6.0*pn[bin-6];
            dpn[bin] /= dR;
         }
         else
         {
            dpn[bin] = -1.0/60.0*pn[bin-3] + 3.0/20.0*pn[bin-2] - 3.0/4.0*pn[bin-1]
                       + 3.0/4.0*pn[bin+1] - 3.0/20.0*pn[bin+2] + 1.0/60.0*pn[bin+3];
            dpn[bin] /= dR;
         }

         integr += pow(r[bin],3.0) * dpn[bin] * dR;
      }

      return integr;
   }


   typedef struct {
      double dl,dg,re,r0,pg_cut,work_cut,eff_st_cut,pg_full,work_full,eff_st_full;
   } Vec_dens_st;

   Vec_dens_st dens_to_st ( double dR, int maxBin, double dl, double dg, double r0, double ksi,
                            double r[], double dens[], double pn_full[], double pn_cut[] )
   {
      Vec_dens_st dens_st;

      int bin, b_beta;
      double re, r_beta;
      double integr, pg_full, pg_cut;
      double work_full, work_cut, eff_st_full, eff_st_cut;

      integr = integr_r3 ( dR, maxBin, r, dens );

      re = integr / ( dg - dl );
      re = pow(re,0.3333333333); /* nm */
/*
 *    read pressure_N and compute effective surface tension
 *    (truncated and full LJ interaction)
 */
      pg_full = 0.0;
      pg_cut = 0.0;
      r_beta = re + ksi*4.0;
      b_beta = (int)( r_beta / dR );
      for ( bin=b_beta; bin<maxBin; bin++ )
      {
         pg_full += pn_full[bin];
         pg_cut += pn_cut[bin];
      }
      pg_full /= (double)( maxBin - b_beta );
      pg_cut /= (double)( maxBin - b_beta );

      work_full = 0.0;
      work_cut = 0.0;
      for ( bin=0; bin<maxBin; bin++ )
      {
         work_full += ( ( pn_full[bin] - pg_full ) * r[bin] * r[bin] ) * dR;
         work_cut += ( ( pn_cut[bin] - pg_cut ) * r[bin] * r[bin] ) * dR;
      }

      work_full *= pi * 2.0 ;  /* kJ mol^-1 */
      eff_st_full = work_full * 3.0 / ( 4.0 * pi * re * re ) ; /* kJ mol^-1 nm^-2 */

      work_cut *= pi * 2.0 ;  /* kJ mol^-1 */
      eff_st_cut = work_cut * 3.0 / ( 4.0 * pi * re * re ) ; /* kJ mol^-1 nm^-2 */
/*
 *    save and return final results
 */
      dens_st.dl          = dl         ;
      dens_st.dg          = dg         ;
      dens_st.re          = re         ;
      dens_st.r0          = r0         ;
      dens_st.pg_cut      = pg_cut     ;
      dens_st.work_cut    = work_cut   ;
      dens_st.eff_st_cut  = eff_st_cut ;
      dens_st.pg_full     = pg_full    ;
      dens_st.work_full   = work_full  ;
      dens_st.eff_st_full = eff_st_full;

      return dens_st;
   }


   typedef struct {
      double pl,pg,rs,r0,st,wk;
   } Vec_pres_st;

   Vec_pres_st pres_to_st ( double dR, int maxBin, double pl, double pg, double r0, double ksi,
                            double r[], double pn[] )
   {
      Vec_pres_st pres_st;

      int bin;
      double integr;
      double st, rs, wk;

      integr = integr_r3 ( dR, maxBin, r, pn );

      st = integr * ( -1.0 / 8.0 ) * ( pl - pg ) * ( pl - pg );
      st = pow(st,0.3333333333); /* kJ mol^-1 nm^-2 */

      rs = st * 2.0 / ( pl - pg ); /* nm, Young-Laplace equation */

      wk = st * rs * rs * pi * 4.0 / 3.0; /* kJ mol^-1 */

      pres_st.pl = pl ;
      pres_st.pg = pg ;
      pres_st.rs = rs ;
      pres_st.r0 = r0 ;
      pres_st.st = st ;
      pres_st.wk = wk ;

      return pres_st;
   }


