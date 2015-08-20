/*****************************************************************************
 This file is part of the SurfTensCalc program.                                      
 Copyright (C) 2012 Xin Li
                                                                           
 Filename:  intramol.h
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

   const double pi   =   3.14159265358979323846 ;

   typedef struct {
      double x, y, z;
   } Vec_R;

   typedef struct {
      Vec_R LJ, LJc, QQ;
   } Vec_3_LJ_QQ;

   typedef struct {
      Vec_R fi, fj;
   } Vec_2;

   typedef struct {
      Vec_R fi, fj, fk;
   } Vec_3;

   typedef struct {
      Vec_R fi, fj, fk, fl;
   } Vec_4;

/*
 * PBC with respect to j
 */
   Vec_R compute_PBC ( double rxi, double ryi, double rzi,
                       double rxj, double ryj, double rzj,
                       double xbox, double ybox, double zbox )
   {
      Vec_R  PBC;
      PBC.x = 0.0; PBC.y = 0.0; PBC.z = 0.0;
      if      ( rxi-rxj < -0.5*xbox ) { PBC.x =  xbox; }
      else if ( rxi-rxj >  0.5*xbox ) { PBC.x = -xbox; }
      if      ( ryi-ryj < -0.5*ybox ) { PBC.y =  ybox; }
      else if ( ryi-ryj >  0.5*ybox ) { PBC.y = -ybox; }
      if      ( rzi-rzj < -0.5*zbox ) { PBC.z =  zbox; }
      else if ( rzi-rzj >  0.5*zbox ) { PBC.z = -zbox; }
      return PBC;
   }

   Vec_3_LJ_QQ compute_LJ_QQ ( double rxi, double ryi, double rzi, 
                               double rxj, double ryj, double rzj,
                               double fudgeLJ, double fudgeQQ, 
                               double c6, double c12, double rCut,
                               double fQQ, double qi, double qj )
   {
      double rij2, rij;
      double rxij, ryij, rzij;
      double rij6, rij12, fij;
      Vec_3_LJ_QQ  LJ_QQ;
/*
 *    distance_ij
 */
      rxij = rxj - rxi;
      ryij = ryj - ryi;
      rzij = rzj - rzi;
      rij2 = rxij*rxij + ryij*ryij + rzij*rzij;
      rij  = sqrt( rij2 );
/*
 *    LJ force:
 *    U = C12 / R^12 - C6 / R^6
 *    F = 12 * C12 / R^13 - 6 * C6 / R^7
 *    vecF = ( 12 * C12 / R^14 - 6 * C6 / R^8 ) * vecR
 *         = ( 2 * C12 / R^12 - C6 / R^6 ) * 6 / R^2 * vecR
 */
      if ( c6!=0.0 && c12!=0.0 )
      {
         rij6 = rij2 * rij2 * rij2 ;
         rij12 = rij6 * rij6;
         fij = ( c12*2.0/rij12 - c6/rij6 ) * 6.0/rij2 * fudgeLJ;
         LJ_QQ.LJ.x = fij * rxij;
         LJ_QQ.LJ.y = fij * ryij;
         LJ_QQ.LJ.z = fij * rzij;
         if ( rij <= rCut )
         {
            LJ_QQ.LJc.x = fij * rxij;
            LJ_QQ.LJc.y = fij * ryij;
            LJ_QQ.LJc.z = fij * rzij;
         }
         else
         {
            LJ_QQ.LJc.x = 0.0 ;
            LJ_QQ.LJc.y = 0.0 ;
            LJ_QQ.LJc.z = 0.0 ;
         }
      }
      else
      {
         LJ_QQ.LJ.x = 0.0 ;
         LJ_QQ.LJ.y = 0.0 ;
         LJ_QQ.LJ.z = 0.0 ;
         LJ_QQ.LJc.x = 0.0 ;
         LJ_QQ.LJc.y = 0.0 ;
         LJ_QQ.LJc.z = 0.0 ;
      }
/*
 *    electrostatic force
 */
      if ( qi!=0.0 && qj!=0.0 )
      {
         fij  = fQQ * qi * qj / ( rij2 * rij ) * fudgeQQ;
         LJ_QQ.QQ.x = fij * rxij;
         LJ_QQ.QQ.y = fij * ryij;
         LJ_QQ.QQ.z = fij * rzij;
      }
      else
      {
         LJ_QQ.QQ.x = 0.0 ;
         LJ_QQ.QQ.y = 0.0 ;
         LJ_QQ.QQ.z = 0.0 ;
      }

      return LJ_QQ;
   }

   Vec_2 compute_bond_1 ( double rxi, double ryi, double rzi, 
                          double rxj, double ryj, double rzj,
                          double b0, double kb )
   {
      double rxij, ryij, rzij, rij2, rij;
      double fij;
      Vec_2 bond_force;

      rxij = rxj - rxi;
      ryij = ryj - ryi;
      rzij = rzj - rzi;
      rij2 = rxij*rxij + ryij*ryij + rzij*rzij;
      rij  = sqrt( rij2 );
/*
 *    Eij = 1/2 kb (rij-r0)^2
 *    fij = -dEij/drij = -kb(rij-r0) = kb(r0-rij)
 *    Vector_fij = fij * Vector_rij / rij
 */
      fij  = kb * ( b0 - rij ) / rij ;
      bond_force.fj.x = fij * rxij;
      bond_force.fj.y = fij * ryij;
      bond_force.fj.z = fij * rzij;
      bond_force.fi.x = -bond_force.fj.x;
      bond_force.fi.y = -bond_force.fj.y;
      bond_force.fi.z = -bond_force.fj.z;

      return bond_force;
   }

   Vec_2 compute_bond_3 ( double rxi, double ryi, double rzi,
                          double rxj, double ryj, double rzj,
                          double b0, double D, double beta )
   {
      double rxij, ryij, rzij, rij2, rij;
      double fij;
      Vec_2 bond_force;

      rxij = rxj - rxi;
      ryij = ryj - ryi;
      rzij = rzj - rzi;
      rij2 = rxij*rxij + ryij*ryij + rzij*rzij;
      rij  = sqrt( rij2 );
/*
 *    Eij = D*[1.0-exp(-beta*(rij-r0))]^2
 *    fij = -dEij/drij = -D * 2.0*[1.0-exp(-beta*(rij-r0))] *
 *                       exp(-beta*(rij-r0)) * beta
 *    Vector_fij = fij * Vector_rij / rij
 */
      fij = -2.0 * D * beta * exp(-beta*(rij-b0)) *
            ( 1.0 - exp(-beta*(rij-b0)) ) / rij;
      bond_force.fj.x = fij * rxij;
      bond_force.fj.y = fij * ryij;
      bond_force.fj.z = fij * rzij;
      bond_force.fi.x = -bond_force.fj.x;
      bond_force.fi.y = -bond_force.fj.y;
      bond_force.fi.z = -bond_force.fj.z;

      return bond_force;
   }


   double min ( double x, double y )
   {
      double min_xy;
      min_xy = x<y ? x : y ;
      return min_xy;
   }

   double max ( double x, double y )
   {
      double max_xy;
      max_xy = x>y ? x : y ;
      return max_xy;
   }

   Vec_R scalar_x_vector ( double ss, Vec_R vv )
   {
      Vec_R sv;
      sv.x = ss * vv.x;
      sv.y = ss * vv.y;
      sv.z = ss * vv.z;
      return sv;
   }

   Vec_R add_vector ( Vec_R aa, Vec_R bb )
   {
      Vec_R a_b;
      a_b.x = aa.x + bb.x;
      a_b.y = aa.y + bb.y;
      a_b.z = aa.z + bb.z;
      return a_b;
   }

   Vec_R substract_vector ( Vec_R aa, Vec_R bb )
   {
      Vec_R a_b;
      a_b.x = aa.x - bb.x;
      a_b.y = aa.y - bb.y;
      a_b.z = aa.z - bb.z;
      return a_b;
   }

   double inner_product ( Vec_R ri, Vec_R rj )
   {
      double inner_p;
      inner_p = ri.x*rj.x + ri.y*rj.y + ri.z*rj.z;
      return inner_p;
   }

   Vec_R outer_product ( Vec_R ri, Vec_R rj )
   {
      Vec_R outer_p;
/*
 *    see http://en.wikipedia.org/wiki/Cross_product
 */
      outer_p.x = ri.y*rj.z - ri.z*rj.y;
      outer_p.y = ri.z*rj.x - ri.x*rj.z;
      outer_p.z = ri.x*rj.y - ri.y*rj.x;
      return outer_p;
   }

   double dist_2 ( Vec_R r )
   {
      double d2;
      d2 = r.x*r.x + r.y*r.y + r.z*r.z;
      return d2;
   }

   Vec_3 compute_angle_1 ( double rxi, double ryi, double rzi,
                           double rxj, double ryj, double rzj,
                           double rxk, double ryk, double rzk,
                           double a0, double ka )
   {
      Vec_R r21, r32, r21_unit, r32_unit, r23, r23_unit;
      double r21_d2, r21_d, r32_d2, r32_d, r23_d2, r23_d;
      double cos123, theta123, sin123, fij;
      Vec_3 angle_force;

      r21.x = rxi - rxj;
      r21.y = ryi - ryj;
      r21.z = rzi - rzj;
      r21_d2 = dist_2( r21 );
      r21_d  = sqrt( r21_d2 );

      r21_unit.x = r21.x / r21_d;
      r21_unit.y = r21.y / r21_d;
      r21_unit.z = r21.z / r21_d;

      r32.x = rxj - rxk;
      r32.y = ryj - ryk;
      r32.z = rzj - rzk;
      r32_d2 = dist_2( r32 );
      r32_d  = sqrt( r32_d2 );

      r32_unit.x = r32.x / r32_d;
      r32_unit.y = r32.y / r32_d;
      r32_unit.z = r32.z / r32_d;

      r23.x = -r32.x;
      r23.y = -r32.y;
      r23.z = -r32.z;
      r23_d2 = r32_d2;
      r23_d  = r32_d;

      cos123 = inner_product(r21,r23) / ( r21_d * r23_d );
      theta123 = acos( cos123 );
      sin123 = sin( theta123 );
      a0 *= pi/180.0;

      fij = ka * ( a0 - theta123 );
      angle_force.fi.x = ( cos123 * r21_unit.x + r32_unit.x ) / ( sin123 * r21_d ) * fij;
      angle_force.fi.y = ( cos123 * r21_unit.y + r32_unit.y ) / ( sin123 * r21_d ) * fij;
      angle_force.fi.z = ( cos123 * r21_unit.z + r32_unit.z ) / ( sin123 * r21_d ) * fij;
      angle_force.fk.x = -( cos123 * r32_unit.x + r21_unit.x ) / ( sin123 * r32_d ) * fij;
      angle_force.fk.y = -( cos123 * r32_unit.y + r21_unit.y ) / ( sin123 * r32_d ) * fij;
      angle_force.fk.z = -( cos123 * r32_unit.z + r21_unit.z ) / ( sin123 * r32_d ) * fij;
      angle_force.fj.x = -angle_force.fi.x - angle_force.fk.x;
      angle_force.fj.y = -angle_force.fi.y - angle_force.fk.y;
      angle_force.fj.z = -angle_force.fi.z - angle_force.fk.z;

      return angle_force;
   }

   Vec_4 compute_dihedral_3 ( double rxi, double ryi, double rzi,
                              double rxj, double ryj, double rzj,
                              double rxk, double ryk, double rzk,
                              double rxl, double ryl, double rzl,
                              double c0, double c1, double c2, 
                              double c3, double c4, double c5 )
   {
      Vec_R  r21, r32, r43, r21_unit, r32_unit, r43_unit;
      double r21_d2, r21_d, r32_d2, r32_d, r43_d2, r43_d;
      Vec_R  r23, r34;
      double r23_d2, r23_d, r34_d2, r34_d;
      double cos123, theta123, sin123;
      double cos234, theta234, sin234;
      Vec_R  r_21_32, r_43_32;
      double const_1, const_4;
      double cos_phi, c123, b432, fij;
      Vec_4  dihedral_force;

      r21.x = rxi - rxj;
      r21.y = ryi - ryj;
      r21.z = rzi - rzj;
      r21_d2 = dist_2( r21 );
      r21_d  = sqrt( r21_d2 );

      r21_unit.x = r21.x / r21_d;
      r21_unit.y = r21.y / r21_d;
      r21_unit.z = r21.z / r21_d;

      r32.x = rxj - rxk;
      r32.y = ryj - ryk;
      r32.z = rzj - rzk;
      r32_d2 = dist_2( r32 );
      r32_d  = sqrt( r32_d2 );

      r32_unit.x = r32.x / r32_d;
      r32_unit.y = r32.y / r32_d;
      r32_unit.z = r32.z / r32_d;

      r43.x = rxk - rxl;
      r43.y = ryk - ryl;
      r43.z = rzk - rzl;
      r43_d2 = dist_2( r43 );
      r43_d  = sqrt( r43_d2 );

      r43_unit.x = r43.x / r43_d;
      r43_unit.y = r43.y / r43_d;
      r43_unit.z = r43.z / r43_d;

      r23.x = -r32.x;
      r23.y = -r32.y;
      r23.z = -r32.z;
      r23_d2 = r32_d2;
      r23_d  = r32_d;

      r34.x = -r43.x;
      r34.y = -r43.y;
      r34.z = -r43.z;
      r34_d2 = r43_d2;
      r34_d  = r43_d;

      cos123 = inner_product(r21,r23) / ( r21_d * r23_d );
      theta123 = acos( cos123 );
      sin123 = sin( theta123 );

      cos234 = inner_product(r32,r34) / ( r32_d * r34_d );
      theta234 = acos( cos234 );
      sin234 = sin( theta234 );

      cos_phi = ( cos123 * cos234 - inner_product(r21_unit,r43_unit) )
                / ( sin123 * sin234 );
/*
 *    Note: in RB function the polymer convention is used
 *    psi = phi - 180(degree)
 *    coefficient (-1)^n is multiplied with c_n
 */
      fij = -( -c1 + c2 * 2.0 * cos_phi 
                   - c3 * 3.0 * pow(cos_phi,2.0) 
                   + c4 * 4.0 * pow(cos_phi,3.0) 
                   - c5 * 5.0 * pow(cos_phi,4.0) );

      r_21_32 = outer_product(r21_unit,r32_unit);
      r_43_32 = outer_product(r43_unit,r32_unit);

      const_1 = inner_product(r43_unit,r_21_32)
                / ( r21_d * pow(sin123,3.0) * sin234);
      const_4 = inner_product(r43_unit,r_21_32)
                / ( r43_d * pow(sin234,3.0) * sin123 );

      c123 = r21_d * cos123 / r32_d - 1.0;
      b432 = r43_d * cos234 / r32_d;

      dihedral_force.fi.x = -const_1 * r_21_32.x * fij;
      dihedral_force.fi.y = -const_1 * r_21_32.y * fij;
      dihedral_force.fi.z = -const_1 * r_21_32.z * fij;

      dihedral_force.fl.x = -const_4 * r_43_32.x * fij;
      dihedral_force.fl.y = -const_4 * r_43_32.y * fij;
      dihedral_force.fl.z = -const_4 * r_43_32.z * fij;

      dihedral_force.fj.x = c123 * dihedral_force.fi.x
                          - b432 * dihedral_force.fl.x;
      dihedral_force.fj.y = c123 * dihedral_force.fi.y
                          - b432 * dihedral_force.fl.y;
      dihedral_force.fj.z = c123 * dihedral_force.fi.z
                          - b432 * dihedral_force.fl.z;

      dihedral_force.fk.x = -dihedral_force.fi.x
                           - dihedral_force.fj.x
                           - dihedral_force.fl.x;
      dihedral_force.fk.y = -dihedral_force.fi.y
                           - dihedral_force.fj.y
                           - dihedral_force.fl.y;
      dihedral_force.fk.z = -dihedral_force.fi.z
                           - dihedral_force.fj.z
                           - dihedral_force.fl.z;

      return dihedral_force;
   }

   Vec_4 compute_dihedral_9 ( double rxi, double ryi, double rzi,
                              double rxj, double ryj, double rzj,
                              double rxk, double ryk, double rzk,
                              double rxl, double ryl, double rzl,
                              double phi0, double kphi, int n ) 
   {
      double c0, c1, c2, c3, c4, c5;
      Vec_4 dihedral_force;

      c0 = 0.0; c1 = 0.0; c2 = 0.0;
      c3 = 0.0; c4 = 0.0; c5 = 0.0;

      if ( n==1 )
      {
         if ( phi0==0.0 )
         {
            c0 = kphi;
            c1 = kphi;
         }
         else if ( phi0==180.0 )
         {
            c0 = kphi;
            c1 = kphi * (-1.0);
         }
         else
         {
            printf ("<> Error: phi0 is not equal to 0 or pi!\n");
            printf ("   phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
            exit(1);
         }
      }
      else if ( n==2 )
      {
         if ( phi0==0.0 )
         {
            c2 = kphi * 2.0;
         }
         else if ( phi0==180.0 )
         {
            c0 = kphi * 2.0;
            c2 = kphi * (-2.0);
         }
         else
         {
            printf ("<> Error: phi0 is not equal to 0 or pi!\n");
            printf ("   phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
            exit(1);
         }
      }
      else if ( n==3 )
      {
         if ( phi0==0.0 )
         {
            c0 = kphi;
            c1 = kphi * (-3.0);
            c3 = kphi * 4.0;
         }
         else if ( phi0==180.0 )
         {
            c0 = kphi;
            c1 = kphi * 3.0;
            c3 = kphi * (-4.0);
         }
         else
         {
            printf ("<> Error: phi0 is not equal to 0 or pi!\n");
            printf ("   phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
            exit(1);
         }
      }
      else if ( n==4 )
      {
         if ( phi0==0.0 )
         {
            c0 = kphi * 2.0;
            c2 = kphi * (-8.0);
            c4 = kphi * 8.0;
         }
         else if ( phi0==180.0 )
         {
            c2 = kphi * 8.0;
            c4 = kphi * (-8.0);
         }
         else
         {
            printf ("<> Error: phi0 is not equal to 0 or pi!\n");
            printf ("   phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
            exit(1);
         }
      }
      else
      {
         printf ("<> Error: n is not equal to 1, 2, 3 or 4!\n");
         printf ("   phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
         exit(1);
      }

      c1 *= -1.0;
      c3 *= -1.0;
      c5 *= -1.0;

      dihedral_force = compute_dihedral_3 ( rxi, ryi, rzi,
                                            rxj, ryj, rzj,
                                            rxk, ryk, rzk,
                                            rxl, ryl, rzl,
                                            c0,  c1,  c2,
                                            c3,  c4,  c5 );

      return dihedral_force;
   }









   void compute_pres_BD_1 ( int maxBin, double dR, double InvD,
                            double inv_dV[], double thrshd,
                            Vec_R vec_ri, double rci, 
                            Vec_R vec_rj, double rcj,
                            Vec_R vec_fi )
   {
      Vec_R vec_rji, vec_rji_rj;
      double cij2, cij, dR2rji2, s, omega2, omega;
      int N_i, N_j, N_0, sFlag, NN, nn, nSlice, bin;
      double *L, *Sigma;
      double ffij_type_1, incr_ff_N1, incr_ff_T1;

      vec_rji.x = vec_ri.x - vec_rj.x;
      vec_rji.y = vec_ri.y - vec_rj.y;
      vec_rji.z = vec_ri.z - vec_rj.z;
      cij2 = dist_2( vec_rji );
      cij = sqrt( cij2 );

      ffij_type_1 = inner_product( vec_fi, vec_rji ) / cij2;

      dR2rji2 = dR*dR*cij2;
      s = inner_product( vec_rji, vec_rj );
      vec_rji_rj = outer_product( vec_rji, vec_rj );
      omega2 = dist_2( vec_rji_rj );
      omega = sqrt( omega2 );

      N_i = (int)(rci*InvD);
      N_j = (int)(rcj*InvD);
      N_0 = (int)( omega/cij * InvD );
      if ( N_0 <= maxBin-1 )
      {
         if ( s < -cij2 )
         {
            sFlag = -1;
            NN = N_j - N_i;
         }
         else if( s <= 0.0 )
         {
            sFlag = 0;
            NN = N_i + N_j - N_0*2;
         }
         else
         {
            sFlag = 1;
            NN = N_i - N_j;
         }
 
         if ( NN < 0 )
         {
            if ( sFlag == 0 )
            {
               printf ("<> Error: NN<0!\n   N_i = %d, N_j = %d, N_0 = %d, NN = %d, sFlag = %d\n",
                        N_i, N_j, N_0, NN, sFlag );
            }
            else
            {
               printf ("<> Error: NN<0!\n   N_i = %d, N_j = %d, NN = %d, sFlag = %d\n",
                        N_i, N_j, NN, sFlag );
            }
            exit(1);
         }
 
         L = malloc( (NN+2) * sizeof(double) );
         Sigma = malloc( (NN+2) * sizeof(double) );
 
         L[0] = rcj*InvD ;
         Sigma[0] = s ;
         nn = 1;
 
         if ( sFlag == -1 )
         {
            for( nSlice = N_j; nSlice > N_i; nSlice-- )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=-sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
         else if ( sFlag == 1 )
         {
            for( nSlice = N_j+1; nSlice <= N_i; nSlice++ )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
 
         else
         {
            for( nSlice = N_j; nSlice > N_0; nSlice-- )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=-sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
            for( nSlice = N_0+1; nSlice <= N_i; nSlice++ )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
         L[NN+1] = rci*InvD ;
         Sigma[NN+1] = cij2 + s ;
 
         if ( nn != NN+1 )
         {
            printf ("<> Error: nn!=NN+1\nnn = %d, NN = %d\n",nn,NN);
            printf ("   N_i = %d, N_j = %d\n",N_i,N_j);
            exit(1);
         }
 
         for ( nn=0; nn<NN+1; nn++ )
         {
            if ( L[nn] == L[nn+1] )
            {
               bin = (int)( 0.5*(L[nn]+L[nn+1]) ) - 1 ;
            }
            else
            {
               bin = (int)( 0.5*(L[nn]+L[nn+1]) );
            }

            if ( bin >= maxBin )
            {
               continue;
            }
 
            if ( omega < thrshd )
            {
               incr_ff_N1 = ( Sigma[nn+1]
                              - Sigma[nn] ) * inv_dV[bin];
               incr_ff_T1 = 0.0;
            }
            else
            {
               incr_ff_N1 = ( Sigma[nn+1] - omega*atan2(Sigma[nn+1],omega)
                              - Sigma[nn] + omega*atan2(Sigma[nn],omega) ) * inv_dV[bin];
               incr_ff_T1 = 0.5 * ( omega*atan2(Sigma[nn+1],omega)
                                    - omega*atan2(Sigma[nn],omega) ) * inv_dV[bin];
            }
 
            pN_U_BD1[bin] += incr_ff_N1 * ffij_type_1 ;

            pT_U_BD1[bin] += incr_ff_T1 * ffij_type_1 ;
         }

         free ( L );
         free ( Sigma );
      }

   }



   void compute_pres_14_1 ( int maxBin, double dR, double InvD,
                            double inv_dV[], double thrshd,
                            Vec_R vec_ri, double rci, 
                            Vec_R vec_rj, double rcj,
                            Vec_R vec_fi )
   {
      Vec_R vec_rji, vec_rji_rj;
      double cij2, cij, dR2rji2, s, omega2, omega;
      int N_i, N_j, N_0, sFlag, NN, nn, nSlice, bin;
      double *L, *Sigma;
      double ffij_type_1, incr_ff_N1, incr_ff_T1;

      vec_rji.x = vec_ri.x - vec_rj.x;
      vec_rji.y = vec_ri.y - vec_rj.y;
      vec_rji.z = vec_ri.z - vec_rj.z;
      cij2 = dist_2( vec_rji );
      cij = sqrt( cij2 );

      ffij_type_1 = inner_product( vec_fi, vec_rji ) / cij2;

      dR2rji2 = dR*dR*cij2;
      s = inner_product( vec_rji, vec_rj );
      vec_rji_rj = outer_product( vec_rji, vec_rj );
      omega2 = dist_2( vec_rji_rj );
      omega = sqrt( omega2 );

      N_i = (int)(rci*InvD);
      N_j = (int)(rcj*InvD);
      N_0 = (int)( omega/cij * InvD );
      if ( N_0 <= maxBin-1 )
      {
         if ( s < -cij2 )
         {
            sFlag = -1;
            NN = N_j - N_i;
         }
         else if( s <= 0.0 )
         {
            sFlag = 0;
            NN = N_i + N_j - N_0*2;
         }
         else
         {
            sFlag = 1;
            NN = N_i - N_j;
         }
 
         if ( NN < 0 )
         {
            if ( sFlag == 0 )
            {
               printf ("<> Error: NN<0!\n   N_i = %d, N_j = %d, N_0 = %d, NN = %d, sFlag = %d\n",
                        N_i, N_j, N_0, NN, sFlag );
            }
            else
            {
               printf ("<> Error: NN<0!\n   N_i = %d, N_j = %d, NN = %d, sFlag = %d\n",
                        N_i, N_j, NN, sFlag );
            }
            exit(1);
         }
 
         L = malloc( (NN+2) * sizeof(double) );
         Sigma = malloc( (NN+2) * sizeof(double) );
 
         L[0] = rcj*InvD ;
         Sigma[0] = s ;
         nn = 1;
 
         if ( sFlag == -1 )
         {
            for( nSlice = N_j; nSlice > N_i; nSlice-- )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=-sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
         else if ( sFlag == 1 )
         {
            for( nSlice = N_j+1; nSlice <= N_i; nSlice++ )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
 
         else
         {
            for( nSlice = N_j; nSlice > N_0; nSlice-- )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=-sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
            for( nSlice = N_0+1; nSlice <= N_i; nSlice++ )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
         L[NN+1] = rci*InvD ;
         Sigma[NN+1] = cij2 + s ;
 
         if ( nn != NN+1 )
         {
            printf ("<> Error: nn!=NN+1\nnn = %d, NN = %d\n",nn,NN);
            printf ("   N_i = %d, N_j = %d\n",N_i,N_j);
            exit(1);
         }
 
         for ( nn=0; nn<NN+1; nn++ )
         {
            if ( L[nn] == L[nn+1] )
            {
               bin = (int)( 0.5*(L[nn]+L[nn+1]) ) - 1 ;
            }
            else
            {
               bin = (int)( 0.5*(L[nn]+L[nn+1]) );
            }

            if ( bin >= maxBin )
            {
               continue;
            }
 
            if ( omega < thrshd )
            {
               incr_ff_N1 = ( Sigma[nn+1]
                              - Sigma[nn] ) * inv_dV[bin];
               incr_ff_T1 = 0.0;
            }
            else
            {
               incr_ff_N1 = ( Sigma[nn+1] - omega*atan2(Sigma[nn+1],omega)
                              - Sigma[nn] + omega*atan2(Sigma[nn],omega) ) * inv_dV[bin];
               incr_ff_T1 = 0.5 * ( omega*atan2(Sigma[nn+1],omega)
                                    - omega*atan2(Sigma[nn],omega) ) * inv_dV[bin];
            }
 
            pN_U_14[bin] += incr_ff_N1 * ffij_type_1 ;

            pT_U_14[bin] += incr_ff_T1 * ffij_type_1 ;
         }

         free ( L );
         free ( Sigma );
      }

   }






   void compute_pres_BD_2 ( int maxBin, double dR, double InvD,
                            double inv_dV[], double thrshd,
                            Vec_R vec_ri, double rci, 
                            Vec_R vec_rj, double rcj,
                            Vec_R vec_fi )
   {
      Vec_R vec_rji, vec_rji_rj;
      double cij2, cij, dR2rji2, s, omega2, omega;
      int N_i, N_j, N_0, sFlag, NN, nn, nSlice, bin;
      double *L;
      double ffij_type_2, incr_ff_N2, incr_ff_T2;

      vec_rji.x = vec_ri.x - vec_rj.x;
      vec_rji.y = vec_ri.y - vec_rj.y;
      vec_rji.z = vec_ri.z - vec_rj.z;
      cij2 = dist_2( vec_rji );
      cij = sqrt( cij2 );

      ffij_type_2 = inner_product( vec_fi, vec_rj );

      dR2rji2 = dR*dR*cij2;
      s = inner_product( vec_rji, vec_rj );
      vec_rji_rj = outer_product( vec_rji, vec_rj );
      omega2 = dist_2( vec_rji_rj );
      omega = sqrt( omega2 );

      N_i = (int)(rci*InvD);
      N_j = (int)(rcj*InvD);
      N_0 = (int)( omega/cij * InvD );
      if ( N_0 <= maxBin-1 )
      {
         if ( s < -cij2 )
         {
            sFlag = -1;
            NN = N_j - N_i;
         }
         else if( s <= 0.0 )
         {
            sFlag = 0;
            NN = N_i + N_j - N_0*2;
         }
         else
         {
            sFlag = 1;
            NN = N_i - N_j;
         }
 
         if ( NN < 0 )
         {
            if ( sFlag == 0 )
            {
               printf ("<> Error: NN<0!\n   N_i = %d, N_j = %d, N_0 = %d, NN = %d, sFlag = %d\n",
                        N_i, N_j, N_0, NN, sFlag );
            }
            else
            {
               printf ("<> Error: NN<0!\n   N_i = %d, N_j = %d, NN = %d, sFlag = %d\n",
                        N_i, N_j, NN, sFlag );
            }
            exit(1);
         }
 
         L = malloc( (NN+2) * sizeof(double) );
 
         L[0] = rcj*InvD ;
         nn = 1;
 
         if ( sFlag == -1 )
         {
            for( nSlice = N_j; nSlice > N_i; nSlice-- )
            {
               L[nn] = (double)nSlice;
               nn++;
            }
         }
         else if ( sFlag == 1 )
         {
            for( nSlice = N_j+1; nSlice <= N_i; nSlice++ )
            {
               L[nn] = (double)nSlice;
               nn++;
            }
         }
 
         else
         {
            for( nSlice = N_j; nSlice > N_0; nSlice-- )
            {
               L[nn] = (double)nSlice;
               nn++;
            }
            for( nSlice = N_0+1; nSlice <= N_i; nSlice++ )
            {
               L[nn] = (double)nSlice;
               nn++;
            }
         }
         L[NN+1] = rci*InvD ;
 
         if ( nn != NN+1 )
         {
            printf ("<> Error: nn!=NN+1\nnn = %d, NN = %d\n",nn,NN);
            printf ("   N_i = %d, N_j = %d\n",N_i,N_j);
            exit(1);
         }
 
         for ( nn=0; nn<NN+1; nn++ )
         {
            if ( L[nn] == L[nn+1] )
            {
               bin = (int)( 0.5*(L[nn]+L[nn+1]) ) - 1 ;
            }
            else
            {
               bin = (int)( 0.5*(L[nn]+L[nn+1]) );
            }

            if ( bin >= maxBin )
            {
               continue;
            }
 
            if ( rci<thrshd || rcj<thrshd )
            {
               incr_ff_N2 = 0.0;
               incr_ff_T2 = 0.0;
            }
            else
            {
               incr_ff_N2 = log( L[nn+1] / L[nn] ) * inv_dV[bin];
               incr_ff_T2 = -0.5 * log( L[nn+1] / L[nn] ) * inv_dV[bin];
            }
 
            pN_U_BD2[bin] += incr_ff_N2 * ffij_type_2 ;

            pT_U_BD2[bin] += incr_ff_T2 * ffij_type_2 ;
         }

         free ( L );
      }

   }




   void compute_pres_NB_1 ( int maxBin, double dR, double InvD,
                            double inv_dV[], double thrshd,
                            Vec_R vec_ri, double rci, 
                            Vec_R vec_rj, double rcj, 
                            Vec_3_LJ_QQ LJ_QQ )
   {
      Vec_R vec_rji, vec_rji_rj;
      double cij2, cij, dR2rji2, s, omega2, omega;
      int N_i, N_j, N_0, sFlag, NN, nn, nSlice, bin;
      double *L, *Sigma;
      double ffij_LJ_1, ffij_LJc_1, ffij_QQ_1, incr_ff_N1, incr_ff_T1;

      vec_rji.x = vec_ri.x - vec_rj.x;
      vec_rji.y = vec_ri.y - vec_rj.y;
      vec_rji.z = vec_ri.z - vec_rj.z;
      cij2 = dist_2( vec_rji );
      cij = sqrt( cij2 );

      ffij_LJ_1  = inner_product( LJ_QQ.LJ,  vec_rji ) / cij2;
      ffij_LJc_1 = inner_product( LJ_QQ.LJc, vec_rji ) / cij2;
      ffij_QQ_1  = inner_product( LJ_QQ.QQ,  vec_rji ) / cij2;

      dR2rji2 = dR*dR*cij2;
      s = inner_product( vec_rji, vec_rj );
      vec_rji_rj = outer_product( vec_rji, vec_rj );
      omega2 = dist_2( vec_rji_rj );
      omega = sqrt( omega2 );

      N_i = (int)(rci*InvD);
      N_j = (int)(rcj*InvD);
      N_0 = (int)( omega/cij * InvD );
      if ( N_0 <= maxBin-1 )
      {
         if ( s < -cij2 )
         {
            sFlag = -1;
            NN = N_j - N_i;
         }
         else if( s <= 0.0 )
         {
            sFlag = 0;
            NN = N_i + N_j - N_0*2;
         }
         else
         {
            sFlag = 1;
            NN = N_i - N_j;
         }
 
         if ( NN < 0 )
         {
            if ( sFlag == 0 )
            {
               printf ("<> Error: NN<0!\n   N_i = %d, N_j = %d, N_0 = %d, NN = %d, sFlag = %d\n",
                        N_i, N_j, N_0, NN, sFlag );
            }
            else
            {
               printf ("<> Error: NN<0!\n   N_i = %d, N_j = %d, NN = %d, sFlag = %d\n",
                        N_i, N_j, NN, sFlag );
            }
            exit(1);
         }
 
         L = malloc( (NN+2) * sizeof(double) );
         Sigma = malloc( (NN+2) * sizeof(double) );
 
         L[0] = rcj*InvD ;
         Sigma[0] = s ;
         nn = 1;
 
         if ( sFlag == -1 )
         {
            for( nSlice = N_j; nSlice > N_i; nSlice-- )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=-sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
         else if ( sFlag == 1 )
         {
            for( nSlice = N_j+1; nSlice <= N_i; nSlice++ )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
 
         else
         {
            for( nSlice = N_j; nSlice > N_0; nSlice-- )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=-sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
            for( nSlice = N_0+1; nSlice <= N_i; nSlice++ )
            {
               L[nn] = (double)nSlice;
               Sigma[nn]=sqrt(dR2rji2*(nSlice*nSlice)-omega2);
               nn++;
            }
         }
         L[NN+1] = rci*InvD ;
         Sigma[NN+1] = cij2 + s ;
 
         if ( nn != NN+1 )
         {
            printf ("<> Error: nn!=NN+1\nnn = %d, NN = %d\n",nn,NN);
            printf ("   N_i = %d, N_j = %d\n",N_i,N_j);
            exit(1);
         }
 
         for ( nn=0; nn<NN+1; nn++ )
         {
            if ( L[nn] == L[nn+1] )
            {
               bin = (int)( 0.5*(L[nn]+L[nn+1]) ) - 1 ;
            }
            else
            {
               bin = (int)( 0.5*(L[nn]+L[nn+1]) );
            }

            if ( bin >= maxBin )
            {
               continue;
            }

/*
if(bin>=maxBin){
printf("bin>=maxBin!\n");
printf("bin=%d maxBin=%d\n",bin,maxBin);
printf("N_i=%d N_0=%d N_j=%d\n",N_i,N_0,N_j);
printf("L[nn]=%f L[nn+1]=%f\n",L[nn],L[nn+1]);
exit(1);
}
*/
 
            if ( omega < thrshd )
            {
               incr_ff_N1 = ( Sigma[nn+1]
                              - Sigma[nn] ) * inv_dV[bin];
               incr_ff_T1 = 0.0;
            }
            else
            {
               incr_ff_N1 = ( Sigma[nn+1] - omega*atan2(Sigma[nn+1],omega)
                              - Sigma[nn] + omega*atan2(Sigma[nn],omega) ) * inv_dV[bin];
               incr_ff_T1 = 0.5 * ( omega*atan2(Sigma[nn+1],omega)
                                    - omega*atan2(Sigma[nn],omega) ) * inv_dV[bin];
            }
 
            pN_U_LJ[bin] += incr_ff_N1 * ffij_LJ_1 ;

            pT_U_LJ[bin] += incr_ff_T1 * ffij_LJ_1 ;

            pN_U_LJc[bin] += incr_ff_N1 * ffij_LJc_1 ;

            pT_U_LJc[bin] += incr_ff_T1 * ffij_LJc_1 ;

            pN_U_QQ[bin] += incr_ff_N1 * ffij_QQ_1 ;

            pT_U_QQ[bin] += incr_ff_T1 * ffij_QQ_1 ;
         }

         free ( L );
         free ( Sigma );
      }

   }


