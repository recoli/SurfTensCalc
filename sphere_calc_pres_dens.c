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
 *    Serial C program to compute surface tension of droplet        *
 *    by Xin Li, TheoChem & Biology, KTH, Stockholm, Sweden         *
 *                      Version 2012-03-16                          *
 ********************************************************************/

   #include <math.h>
   #include <time.h>
   #include <stdio.h>
   #include <stdlib.h>
   #include <xdrfile/xdrfile_trr.h>

   double *pN_U_LJ, *pN_U_LJc, *pN_U_QQ;
   double *pN_U_14, *pN_U_BD1, *pN_U_BD2;

   double *pT_U_LJ, *pT_U_LJc, *pT_U_QQ;
   double *pT_U_14, *pT_U_BD1, *pT_U_BD2;

   time_t start_t, end_t;
   int delta_time;

   #include "intramol.h"

/********************************************************************
 *    Main program                                                  *
 ********************************************************************/

   int main( int argc, char *argv[] )
   {
/*
 *    There should be 5 arguments: nFrames, nStart, temperature, rCut, boxSize
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

      start_t = time(NULL);
      printf ( ">>> sphere_calc_pres_dens.x  version 1.8.1 \n" );
      printf ( ">>> Calculation started at %s", ctime(&start_t) );
      printf ( ">>> %s %s %s %s %s %s\n",
               argv[0], argv[1], argv[2], argv[3], argv[4], argv[5] );


/********************************************************************
 *    Define variables                                              *
 ********************************************************************/

      Vec_3_LJ_QQ   LJ_QQ;

      Vec_2   bond_force;
      double  b0, kb, D, beta;

      Vec_3   angle_force;
      double  a0, ka;

      Vec_4   dihedral_force;
      int    k, l; 
      double c0, c1, c2, c3, c4, c5; 
      double phi0, kphi; 
      int    n;


      typedef struct {
         double charge, mass, sigma, epsilon;
         int atomtype, residue;
      } AtomParam;

      AtomParam   **atom_param;

/*
 *    define structure types
 */
      typedef struct {
         int atom_i, atom_j; 
         int res_i, res_j; 
         double b0, kb; 
      } Bond_1;

      typedef struct {
         int atom_i, atom_j;
         int res_i, res_j;
         double b0, D, beta;
      } Bond_3;

      Bond_1   **bond_1;
      Bond_3   **bond_3;
      int      iBond, *nBonds;
      int      **bond_funct;

      typedef struct {
         int atom_i, atom_j;
         int res_i, res_j;
      } Pair_1;

      typedef struct {
         int atom_i, atom_j;
         int res_i, res_j;
         double fudgeQQ, q_i, q_j, V, W;
      } Pair_2;

      Pair_1   **pair_1;
      Pair_2   **pair_2;
      int      iPair, *nPairs;
      int      **pair_funct;
      double   pair_fudgeQQ;

      typedef struct {
         int atom_i, atom_j, atom_k;
         int res_i, res_j, res_k;
         double a0, ka;
      } Angle_1;

      Angle_1   **angle_1;
      int       iAngle, *nAngles;
      int       **angle_funct;

      typedef struct {
         int atom_i, atom_j, atom_k, atom_l;
         int res_i, res_j, res_k, res_l;
         double c0, c1, c2, c3, c4, c5;
      } Dihedral_3;

      typedef struct {
         int atom_i, atom_j, atom_k, atom_l;
         int res_i, res_j, res_k, res_l;
         double phi0, kphi;
         int n;
      } Dihedral_9;

      Dihedral_3   **dihedral_3;
      Dihedral_9   **dihedral_9;
      int   iDihedral, *nDihedrals;
      int   **dihedral_funct;

      int atom_i, atom_j, atom_k, atom_l;
      int res_i, res_j, res_k, res_l;
      int funct;

      int ***exclude;
/*
 *    number of frames in trajectory file 
 */
      int   frame, nFrames, nStart;
      nFrames = atoi(argv[1]);
      nStart = atoi(argv[2]);

      if ( nFrames<=nStart ) 
      {
         printf( "Error: nFrames is no larger than nStart!\n" ) ;
         printf( "nFrames = %d, nStart = %d\n", nFrames, nStart );
         exit(1);
      }

      double temperature;
      temperature = atof(argv[3]); /* K */

      double rCut;
      rCut = atof(argv[4]); /* nm */

      double boxSize;
      boxSize = atof(argv[5]); /* nm */

/*
 *    Constants    
 *
 *    Source: 2010 CODATA
 *    c (speed of light in vacuum):  299792458 m s^-1
 *    mu_0 (magnetic constant):      4*pi * 10^-7 N A^-2 (A=Ampere)
 *    e (elementary charge):         1.602176565 * 10^-19 C 
 *    nA (Avogadro constant):        6.02214129 * 10^23 mol^-1
 *    kB (Boltzmann constant):       1.3806488 * 10^-23 J K^-1
 *
 *    Coulomb constant k=1/(4*pi*epsilon_0)=c^2*mu_0/(4*pi)
 *    fQQ = 2.99792458^2 * 10^16 * 10^-7 N m^2 C^-2
 *        = 2.99792458^2 * 6.02214129 * 1.602176565^2 kJ mol^-1 nm e^-2
 *        = 138.93545782935
 *
 *    PI to 20 decimal places:       3.14159265358979323846
 */
      const double kB   =   1.3806488e-3 ; /* 10^-23 kJ K^-1 */
      const double fQQ  = 138.93545782935 ;    /* kJ mol^-1 nm e^-2 */
      const double pi   =   3.14159265358979323846 ;
      const double nA   =   6.02214129 ;  /* 10^23 mol^-1 */

/*
 *    kT = kB*nA*temperature, in kJ mol^-1
 */
      double kT;
      kT = kB*nA*temperature;

      int comb_rule;
      double fudgeLJ, fudgeQQ;

      int nTypes;
      double **LJ_C6, **LJ_C12;
      int iType, jType;
      double c6, c12;

/*
 *    work of formation
 */
      double work_N_full, work_N_cut, work_T_full, work_T_cut;
/*
 *    coordinates of atoms and COMs, COM = center of mass 
 */
      double *rx, *ry, *rz;
      double *vx, *vy, *vz;
      double *fx, *fy, *fz;
/*
      double **cx, **cy, **cz;
*/
/*
 *    thickness of sperical layer, 0.3 Angstrom = 0.1 Sigma_Oxygen
 */  
      const double dR = 0.05 ;
      double InvD = 1.0/dR;

      int maxBin ;
      maxBin = (int)( ( boxSize*0.5 - 1.0 ) * InvD ) ;


      const double thrshd=1.0e-04;
      Vec_R vec_box, vec_ri, vec_vi, vec_fi;
      Vec_R vec_rj, vec_rji, vec_rji_rj ;

      double vci2, kE_N, kE_T;
      int iia;


      int    bin ;
      double *inv_dV ;
      double *dens ;

      double *pN_K, *pN, *pNc;
      double *pT_K, *pT, *pTc;
/*
 *    densities for different molecules 
 */
      //double **molDens ;
      double ***resDens ;
/*
 *    increment, temporary variable 
 */
      double incr_dens; 
/*
 *    pair pointers for loop   
 */
      long int pair, totalPairs;
      int    *molI, *molJ;
/*
 *    center of mass   
 */
      double xcom, ycom, zcom, msum ;
/*
 *    variables for pair force calculation   
 */
      int    i, j, im, jm, iMin, iMax, jMin, jMax ;
      double mi, qi, qj ;
      double r ;

      double *r_ci;

      double rci, rcj ;


      int x_bin, y_bin, z_bin;
      int x_maxBin, y_maxBin, z_maxBin;
      int *x_dens, *y_dens, *z_dens;
      double x_ref, y_ref, z_ref;
      double x_o, y_o, z_o;

      x_maxBin = (int)( boxSize*2.0 * InvD ) ;
      y_maxBin = (int)( boxSize*2.0 * InvD ) ;
      z_maxBin = (int)( boxSize*2.0 * InvD ) ;

      x_dens = malloc(x_maxBin * sizeof(int));
      y_dens = malloc(y_maxBin * sizeof(int));
      z_dens = malloc(z_maxBin * sizeof(int));


/*
 *    variables for xdrfile_trr  
 */
      int     gmxStep,read_return,gmxNAtoms;
      float   gmxTime,gmxLambda;
      matrix  gmxBox;
      rvec    *x, *v, *f;
      XDRFILE *trr;


/********************************************************************
 *    Reading parameters from param.txt                             *
 ********************************************************************/

      FILE    *file_par, *file_dens, *file_pres;
/*
 *    atom = index of atoms in each molecule   
 *    mol  = index of molecules in the system 
 */
      int     mol, atom, molTypes;
      int     res, ires, jres, km;
/*
 *    atom_in_res = number of atoms in each type of molecule & residue
 *    atom_in_mol = number of atoms in each type of molecule
 *    resNr  = number of residues in each type of molecule
 *    molNr  = number of each type of molecule in the system 
 */
      int     **atom_in_res, *atom_in_mol;
      int     *resNr, *molNr, *iAtom, *iMol, **mini, **maxi;
      char    line[256], tmp[256], filenameTRR[256];
      int     nAtoms, nMols;
      double  *molMass;
      double  **resMass;
/*
 *    open param.txt 
 */
      file_par = fopen("param.txt","r") ;
      if ( NULL==file_par ) {
         printf( "Error: Cannot open file param.txt !\n" ) ;
         exit(1);
      }
      mol = 0;
/*
 *    read comb_rule and fudgeLJ, fudgeQQ (first two lines)
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%d", tmp, &comb_rule );
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%lf%lf", tmp, &fudgeLJ, &fudgeQQ );
/*
 *    read name of trr file
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%s", tmp, filenameTRR );
/*
 *    read number of molecules
 *    atom_in_res: number of atoms in this type of molecule and residue
 *    atom_in_mol = number of atoms in each type of molecule
 *    molNr: number of this type of molecule in the system
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%d", tmp, &molTypes );

      atom_in_res = malloc(molTypes * sizeof(int *));
      atom_in_mol = malloc(molTypes * sizeof(int));

      resNr  = malloc(molTypes * sizeof(int));
      molNr  = malloc(molTypes * sizeof(int));

      for ( mol=0; mol<molTypes; mol++ ) 
      {
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d%d", tmp, &molNr[mol], &resNr[mol] );

         atom_in_res[mol]  = malloc(resNr[mol] * sizeof(int));
         for ( res=0; res<resNr[mol]; res++ )
         {
            atom_in_res[mol][res] = 0;
         }

      }

/*
 *    read atomic information
 */
      atom_param = malloc(molTypes * sizeof(AtomParam *));
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &atom_in_mol[mol] );
         atom_param[mol] = malloc(atom_in_mol[mol] * sizeof(AtomParam));
         for ( atom=0; atom<atom_in_mol[mol]; atom++ ) 
         {
            fgets( line, sizeof( line ), file_par );
            sscanf( line, "%lf%lf%lf%lf%d%d",
                    &(atom_param[mol][atom].charge),
                    &(atom_param[mol][atom].mass),
                    &(atom_param[mol][atom].sigma),
                    &(atom_param[mol][atom].epsilon),
                    &(atom_param[mol][atom].atomtype),
                    &(atom_param[mol][atom].residue) );
            res = atom_param[mol][atom].residue;
            atom_in_res[mol][res]++;
         }
      }
/*
 *    read Lennard-Jones C6 and C12 coefficients
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%d", tmp, &nTypes );
      LJ_C6  = malloc(nTypes * sizeof(double *));
      LJ_C12 = malloc(nTypes * sizeof(double *));
      for ( iType=0; iType<nTypes; iType++ )
      {
         LJ_C6[iType]  = malloc(nTypes * sizeof(double));
         LJ_C12[iType] = malloc(nTypes * sizeof(double));
         for ( jType=iType; jType<nTypes; jType++ )
         {
            fgets( line, sizeof( line ), file_par );
            sscanf( line, "%le%le", &LJ_C6[iType][jType], 
                                    &LJ_C12[iType][jType] );
         }
      }

      for ( iType=0; iType<nTypes-1; iType++ )
      {
         for ( jType=iType+1; jType<nTypes; jType++ )
         {
            LJ_C6[jType][iType] = LJ_C6[iType][jType];
            LJ_C12[jType][iType] = LJ_C12[iType][jType];
         }
      }
/*
 *    read bonded potentials: bond, pair, angle, dihedral
 */
      nBonds = malloc ( molTypes * sizeof( int ) ) ;
      bond_1 = malloc ( molTypes * sizeof( Bond_1 * ) ) ;
      bond_3 = malloc ( molTypes * sizeof( Bond_3 * ) ) ;
      bond_funct = malloc ( molTypes * sizeof( int * ) ) ;

      nPairs = malloc ( molTypes * sizeof( int ) ) ;
      pair_1 = malloc ( molTypes * sizeof( Pair_1 * ) ) ;
      pair_2 = malloc ( molTypes * sizeof( Pair_2 * ) ) ;
      pair_funct = malloc ( molTypes * sizeof( int * ) ) ;

      nAngles = malloc ( molTypes * sizeof( int ) ) ;
      angle_1 = malloc ( molTypes * sizeof( Angle_1 * ) ) ;
      angle_funct = malloc ( molTypes * sizeof( int * ) ) ;

      nDihedrals = malloc ( molTypes * sizeof( int ) ) ;
      dihedral_3 = malloc ( molTypes * sizeof( Dihedral_3 * ) ) ;
      dihedral_9 = malloc ( molTypes * sizeof( Dihedral_9 * ) ) ;
      dihedral_funct = malloc ( molTypes * sizeof( int * ) ) ;

      exclude = malloc ( molTypes * sizeof( int ** ) ) ;

      for ( mol=0; mol<molTypes; mol++ )
      {
/*
 *       exclusions
 */
         exclude[mol] = malloc ( atom_in_mol[mol] * sizeof( int * ) ) ;
         for ( atom_i=0; atom_i<atom_in_mol[mol]; atom_i++ ) 
         {
            exclude[mol][atom_i] = malloc ( atom_in_mol[mol] * sizeof( int ) ) ;
            for ( atom_j=0; atom_j<atom_in_mol[mol]; atom_j++ ) 
            {
               exclude[mol][atom_i][atom_j] = 0;
            }
         }
/*
 *       read bonds
 */
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &nBonds[mol] );

         bond_1[mol] = malloc ( nBonds[mol] * sizeof( Bond_1 ) ) ;
         bond_3[mol] = malloc ( nBonds[mol] * sizeof( Bond_3 ) ) ;
         bond_funct[mol] = malloc ( nBonds[mol] * sizeof( int ) ) ;

         for ( iBond=0; iBond<nBonds[mol]; iBond++ )
         {
            fgets( line, sizeof( line ), file_par );
            sscanf( line, "%d%d%d%d%d",
                    &atom_i, &atom_j, &res_i, &res_j, &funct );
/*
 *          Harmonic bond potential
 */
            if ( funct == 1 )
            {
               sscanf( line, "%d%d%d%d%d%le%le",
                       &(bond_1[mol][iBond].atom_i), &(bond_1[mol][iBond].atom_j),
                       &(bond_1[mol][iBond].res_i), &(bond_1[mol][iBond].res_j),
                       &(bond_funct[mol][iBond]),
                       &(bond_1[mol][iBond].b0), &(bond_1[mol][iBond].kb) );
               exclude[mol][atom_i][atom_j] = 1;
               exclude[mol][atom_j][atom_i] = 1;
            }
/*
 *          Morse bond potential
 */
            else if ( funct == 3 )
            {
               sscanf( line, "%d%d%d%d%d%le%le%le",
                       &(bond_3[mol][iBond].atom_i), &(bond_3[mol][iBond].atom_j),
                       &(bond_3[mol][iBond].res_i), &(bond_3[mol][iBond].res_j),
                       &(bond_funct[mol][iBond]),
                       &(bond_3[mol][iBond].b0), &(bond_3[mol][iBond].D), &(bond_3[mol][iBond].beta) );
               exclude[mol][atom_i][atom_j] = 1;
               exclude[mol][atom_j][atom_i] = 1;
            }
/*
 *          exclusions
 */
            else if ( funct == 5 )
            {
               exclude[mol][atom_i][atom_j] = 1;
               exclude[mol][atom_j][atom_i] = 1;
            }

            else
            {
               printf ( "Error: unsupported bond function type: %d in line:\n%s\n", 
                        funct, line );
               exit(1);
            }
         }
/*
 *       read pairs
 */
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &nPairs[mol] );

         pair_1[mol] = malloc ( nPairs[mol] * sizeof( Pair_1 ) ) ;
         pair_2[mol] = malloc ( nPairs[mol] * sizeof( Pair_2 ) ) ;
         pair_funct[mol] = malloc ( nPairs[mol] * sizeof( int ) ) ;

         for ( iPair=0; iPair<nPairs[mol]; iPair++ )
         {
            fgets( line, sizeof( line ), file_par );
            sscanf( line, "%d%d%d%d%d",
                    &atom_i, &atom_j, &res_i, &res_j, &funct );
            if ( funct == 1 )
            {
               sscanf( line, "%d%d%d%d%d",
                       &(pair_1[mol][iPair].atom_i), &(pair_1[mol][iPair].atom_j),
                       &(pair_1[mol][iPair].res_i), &(pair_1[mol][iPair].res_j),
                       &(pair_funct[mol][iPair]) );
               exclude[mol][atom_i][atom_j] = 1;
               exclude[mol][atom_j][atom_i] = 1;
            }
            else if ( funct == 2 )
            {
               sscanf( line, "%d%d%d%d%d%lf%lf%lf%le%le",
                       &(pair_2[mol][iPair].atom_i), &(pair_2[mol][iPair].atom_j),
                       &(pair_2[mol][iPair].res_i), &(pair_2[mol][iPair].res_j),
                       &(pair_funct[mol][iPair]),
                       &(pair_2[mol][iPair].fudgeQQ),
                       &(pair_2[mol][iPair].q_i), &(pair_2[mol][iPair].q_j),
                       &(pair_2[mol][iPair].V), &(pair_2[mol][iPair].W) );
               exclude[mol][atom_i][atom_j] = 1;
               exclude[mol][atom_j][atom_i] = 1;
            }
            else
            {
               printf ( "Error: unsupported pair function type: %d in line:\n%s\n", 
                        funct, line );
               exit(1);
            }
         }
/*
 *       read angles
 */
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &nAngles[mol] );

         angle_1[mol] = malloc ( nAngles[mol] * sizeof( Angle_1 ) ) ;
         angle_funct[mol] = malloc ( nAngles[mol] * sizeof( int ) ) ;

         for ( iAngle=0; iAngle<nAngles[mol]; iAngle++ )
         {
            fgets( line, sizeof( line ), file_par );
            sscanf( line, "%d%d%d%d%d%d%d",
                    &atom_i, &atom_j, &atom_k, 
                    &res_i, &res_j, &res_k, &funct );
            if ( funct == 1 )
            {
               sscanf( line, "%d%d%d%d%d%d%d%le%le",
                       &(angle_1[mol][iAngle].atom_i), &(angle_1[mol][iAngle].atom_j), 
                       &(angle_1[mol][iAngle].atom_k),
                       &(angle_1[mol][iAngle].res_i), &(angle_1[mol][iAngle].res_j), 
                       &(angle_1[mol][iAngle].res_k),
                       &(angle_funct[mol][iAngle]),
                       &(angle_1[mol][iAngle].a0), &(angle_1[mol][iAngle].ka) );
               exclude[mol][atom_i][atom_j] = 1;
               exclude[mol][atom_j][atom_i] = 1;
               exclude[mol][atom_i][atom_k] = 1;
               exclude[mol][atom_k][atom_i] = 1;
               exclude[mol][atom_j][atom_k] = 1;
               exclude[mol][atom_k][atom_j] = 1;
            }
            else
            {
               printf ( "Error: unsupported angle function type: %d in line:\n%s\n", 
                        funct, line );
               exit(1);
            }
         }
/*
 *       read dihedrals
 */
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &nDihedrals[mol] );

         dihedral_3[mol] = malloc ( nDihedrals[mol] * sizeof( Dihedral_3 ) ) ;
         dihedral_9[mol] = malloc ( nDihedrals[mol] * sizeof( Dihedral_9 ) ) ;
         dihedral_funct[mol] = malloc ( nDihedrals[mol] * sizeof( int ) ) ;

         for ( iDihedral=0; iDihedral<nDihedrals[mol]; iDihedral++ )
         {
            fgets( line, sizeof( line ), file_par );
            sscanf( line, "%d%d%d%d%d%d%d%d%d",
                    &atom_i, &atom_j, &atom_k, &atom_l,
                    &res_i, &res_j, &res_k, &res_l, &funct );
            if ( funct == 3 )
            {
               sscanf( line, "%d%d%d%d%d%d%d%d%d%lf%lf%lf%lf%lf%lf",
                       &(dihedral_3[mol][iDihedral].atom_i), &(dihedral_3[mol][iDihedral].atom_j),
                       &(dihedral_3[mol][iDihedral].atom_k), &(dihedral_3[mol][iDihedral].atom_l),
                       &(dihedral_3[mol][iDihedral].res_i), &(dihedral_3[mol][iDihedral].res_j),
                       &(dihedral_3[mol][iDihedral].res_k), &(dihedral_3[mol][iDihedral].res_l),
                       &(dihedral_funct[mol][iDihedral]),
                       &(dihedral_3[mol][iDihedral].c0), &(dihedral_3[mol][iDihedral].c1),
                       &(dihedral_3[mol][iDihedral].c2), &(dihedral_3[mol][iDihedral].c3),
                       &(dihedral_3[mol][iDihedral].c4), &(dihedral_3[mol][iDihedral].c5) );
               exclude[mol][atom_i][atom_j] = 1;
               exclude[mol][atom_j][atom_i] = 1;
               exclude[mol][atom_i][atom_k] = 1;
               exclude[mol][atom_k][atom_i] = 1;
               exclude[mol][atom_i][atom_l] = 1;
               exclude[mol][atom_l][atom_i] = 1;
               exclude[mol][atom_j][atom_k] = 1;
               exclude[mol][atom_k][atom_j] = 1;
               exclude[mol][atom_j][atom_l] = 1;
               exclude[mol][atom_l][atom_j] = 1;
               exclude[mol][atom_k][atom_l] = 1;
               exclude[mol][atom_l][atom_k] = 1;
            }
            else if ( funct==1 || funct==4 || funct==9 )
            {
               sscanf( line, "%d %d %d %d   %d %d %d %d   %d   %lf %lf %d ",
                       &(dihedral_9[mol][iDihedral].atom_i), &(dihedral_9[mol][iDihedral].atom_j),
                       &(dihedral_9[mol][iDihedral].atom_k), &(dihedral_9[mol][iDihedral].atom_l),
                       &(dihedral_9[mol][iDihedral].res_i), &(dihedral_9[mol][iDihedral].res_j),
                       &(dihedral_9[mol][iDihedral].res_k), &(dihedral_9[mol][iDihedral].res_l),
                       &(dihedral_funct[mol][iDihedral]),
                       &(dihedral_9[mol][iDihedral].phi0), 
                       &(dihedral_9[mol][iDihedral].kphi),
                       &(dihedral_9[mol][iDihedral].n) );
               exclude[mol][atom_i][atom_j] = 1;
               exclude[mol][atom_j][atom_i] = 1;
               exclude[mol][atom_i][atom_k] = 1;
               exclude[mol][atom_k][atom_i] = 1;
               exclude[mol][atom_i][atom_l] = 1;
               exclude[mol][atom_l][atom_i] = 1;
               exclude[mol][atom_j][atom_k] = 1;
               exclude[mol][atom_k][atom_j] = 1;
               exclude[mol][atom_j][atom_l] = 1;
               exclude[mol][atom_l][atom_j] = 1;
               exclude[mol][atom_k][atom_l] = 1;
               exclude[mol][atom_l][atom_k] = 1;
            }
            else
            {
               printf ( "Error: unsupported dihedral function type: %d in line:\n%s\n", 
                        funct, line );
               exit(1);
            }

         }

      }
/*
 *    close param.txt
 */
      fclose( file_par );

/********************************************************************
 *    Check number of molecules and atoms in each molecule          *
 ********************************************************************/

/*
 *    calculate total number of atoms and molecules
 */
      nAtoms = 0;
      nMols = 0;
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         nAtoms += molNr[mol]*atom_in_mol[mol];
         nMols += molNr[mol];
      }
      iAtom = malloc(nAtoms * sizeof(int));
      iMol  = malloc(nAtoms * sizeof(int));
/*
 *    determine mass of each type of molecule
 */
      molMass = malloc(molTypes * sizeof(double));
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         molMass[mol] = 0.0;
         for ( atom=0; atom<atom_in_mol[mol]; atom++ ) 
         {
            molMass[mol] += atom_param[mol][atom].mass;
         }
      }

      resMass = malloc(molTypes * sizeof(double *));
      for ( mol=0; mol<molTypes; mol++ )
      {
         resMass[mol]  = malloc(resNr[mol] * sizeof(double));

         iMin = 0;
         iMax = 0;
         for ( res=0; res<resNr[mol]; res++ )
         {
            resMass[mol][res] = 0.0;
            iMax += atom_in_res[mol][res];
            for ( atom=iMin; atom<iMax; atom++ )
            {
               resMass[mol][res] += atom_param[mol][atom].mass;
            }
            iMin += atom_in_res[mol][res];
         }
      }

/*
 *    assign atom index and mol index for each atom in the system 
 *    i = atom index, im = mol index
 */
      i = 0;
      im = 0;
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         for ( km=0; km<molNr[mol]; km++ ) 
         {
            for ( atom=0; atom<atom_in_mol[mol]; atom++ ) 
            {
               iMol[i] = mol;
               iAtom[i] = atom;
               i++;
            }
            im++;
         }
      }

      if ( nAtoms != i ) 
      {
         printf("Error: Incorrect number of atoms in the system! (1st)\n");
         printf("nAtoms = %d, i = %d\n", nAtoms, i);
         exit(1);
      }

      if ( nMols != im ) 
      {
         printf("Error: Incorrect number of molecules in the system! (1st)\n");
         printf("nMols= %d, im = %d\n", nMols, im);
         exit(1);
      }
/*
 *    assign starting index and ending index for each molecule in the system 
 *    i = atom index, im = mol index
 */

      mini = malloc(nMols * sizeof(int *));
      maxi = malloc(nMols * sizeof(int *));
      im = 0;
      for ( mol=0; mol<molTypes; mol++ )
      {
         for ( km=0; km<molNr[mol]; km++ )
         {
            mini[im] = malloc(resNr[mol] * sizeof(int));
            maxi[im] = malloc(resNr[mol] * sizeof(int));
            im++;
         }
      }

      i = 0;
      im = 0;
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         for ( km=0; km<molNr[mol]; km++ ) 
         {
            for ( res=0; res<resNr[mol]; res++ )
            {
               mini[im][res] = i;
               i += atom_in_res[mol][res];
               maxi[im][res] = i;
            }
            im++;
         }
      }

      if ( nAtoms != i ) 
      {
         printf("Error: Incorrect number of atoms in the system! (2nd)\n");
         printf("nAtoms = %d, i = %d\n", nAtoms, i);
         exit(1);
      }

      if ( nMols != im ) 
      {
         printf("Error: Incorrect number of molecules in the system! (2nd)\n");
         printf("nMols= %d, im = %d\n", nMols, im);
         exit(1);
      }

/********************************************************************
 *    Allocating arrays                                             *
 ********************************************************************/

      rx = malloc( sizeof(double) * nAtoms );
      ry = malloc( sizeof(double) * nAtoms );
      rz = malloc( sizeof(double) * nAtoms );

      vx = malloc( sizeof(double) * nAtoms );
      vy = malloc( sizeof(double) * nAtoms );
      vz = malloc( sizeof(double) * nAtoms );

      fx = malloc( sizeof(double) * nAtoms );
      fy = malloc( sizeof(double) * nAtoms );
      fz = malloc( sizeof(double) * nAtoms );

      r_ci  = malloc( sizeof(double) * nAtoms );
/*
      cx = malloc( sizeof(double *) * nMols );
      cy = malloc( sizeof(double *) * nMols );
      cz = malloc( sizeof(double *) * nMols );
      im = 0;
      for ( mol=0; mol<molTypes; mol++ )
      {
         for ( km=0; km<molNr[mol]; km++ )
         {
            cx[im]  = malloc(resNr[mol] * sizeof(double));
            cy[im]  = malloc(resNr[mol] * sizeof(double));
            cz[im]  = malloc(resNr[mol] * sizeof(double));
            im++;
         }
      }
*/

      inv_dV = malloc( sizeof( double ) * maxBin );

      dens  = malloc( sizeof( double ) * maxBin );

      pN_U_LJ  = malloc( sizeof( double ) * maxBin );
      pN_U_LJc = malloc( sizeof( double ) * maxBin );
      pN_U_QQ  = malloc( sizeof( double ) * maxBin );
      pN_U_14  = malloc( sizeof( double ) * maxBin );
      pN_U_BD1 = malloc( sizeof( double ) * maxBin );
      pN_U_BD2 = malloc( sizeof( double ) * maxBin );
      pN_K = malloc( sizeof( double ) * maxBin );
      pN  = malloc( sizeof( double ) * maxBin );
      pNc = malloc( sizeof( double ) * maxBin );

      pT_U_LJ  = malloc( sizeof( double ) * maxBin );
      pT_U_LJc = malloc( sizeof( double ) * maxBin );
      pT_U_QQ  = malloc( sizeof( double ) * maxBin );
      pT_U_14  = malloc( sizeof( double ) * maxBin );
      pT_U_BD1 = malloc( sizeof( double ) * maxBin );
      pT_U_BD2 = malloc( sizeof( double ) * maxBin );
      pT_K = malloc( sizeof( double ) * maxBin );
      pT  = malloc( sizeof( double ) * maxBin );
      pTc = malloc( sizeof( double ) * maxBin );


      //molDens = malloc(maxBin * sizeof(double *));
      resDens = malloc(maxBin * sizeof(double **));
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         //molDens[bin] = malloc(molTypes * sizeof(double));
         resDens[bin] = malloc(molTypes * sizeof(double *));
         for ( mol=0; mol<molTypes; mol++ )
         {
            resDens[bin][mol] = malloc(resNr[mol] * sizeof(double));
         }
      }

      x = calloc(nAtoms,sizeof(x[0]));
      v = calloc(nAtoms,sizeof(x[0]));
      f = calloc(nAtoms,sizeof(x[0]));

      totalPairs = nMols*(nMols-1)/2 ;
      molI = malloc( sizeof( int ) * totalPairs );
      molJ = malloc( sizeof( int ) * totalPairs );


/********************************************************************
 *    Calculating pressure profile and work of formation            *
 ********************************************************************/

      for ( bin=0; bin<maxBin; bin++ ) 
      {
/*            
 *       volume of the spherical layer 
 */
         inv_dV[bin] = 1.0 / ( pi*4.0/3.0*( pow((double)(bin+1)*dR,3.0) 
                                            - pow((double)bin*dR,3.0) ) );
/*
 *       initialize array: resDens
 */
      }
/*
 *    read trr files for coordinates   
 */
      trr = xdrfile_open (filenameTRR,"r");
      read_trr_natoms (filenameTRR, &gmxNAtoms);
/*
 *    start of loop (frame)   
 */
      for (frame=0; frame<nFrames; frame++) 
      {
/*
 *       read trr file   
 */
         read_return = read_trr (trr,gmxNAtoms,&gmxStep,&gmxTime,&gmxLambda,gmxBox,x,v,f);
         if (read_return!=0) {
            printf ("Error: Failed to read trr file!\n");
            printf ("Frame: %d\n", frame);
            exit(1);
         }
         if (gmxNAtoms!=nAtoms) {
            printf ("Error: Incorrect number of Atoms!\n");
            printf ("gmxNAtoms = %d, nAtoms = %d\n", gmxNAtoms, nAtoms);
            exit(1);
         }
/*
 *       skip if frame < nStart
 */
         if (frame<nStart) 
         { 
            continue; 
         }
/*
 *       assign coordinates
 */
         for ( i=0; i<nAtoms; i++ ) 
         {
            rx[i] = (double)(x[i][0]);
            ry[i] = (double)(x[i][1]);
            rz[i] = (double)(x[i][2]);

            vx[i] = (double)(v[i][0]);
            vy[i] = (double)(v[i][1]);
            vz[i] = (double)(v[i][2]);

            fx[i] = (double)(f[i][0]);
            fy[i] = (double)(f[i][1]);
            fz[i] = (double)(f[i][2]);
         }

         vec_box.x = (double)(gmxBox[0][0]);
         vec_box.y = (double)(gmxBox[1][1]);
         vec_box.z = (double)(gmxBox[2][2]);


/*
 *       find maximal density along three dimensions x,y,z
 *       (corresponding to center of droplet)
 */
         x_maxBin = (int)( vec_box.x * InvD ) ;
         y_maxBin = (int)( vec_box.y * InvD ) ;
         z_maxBin = (int)( vec_box.z * InvD ) ;

         for ( x_bin=0; x_bin<x_maxBin; x_bin++ ) { x_dens[x_bin] = 0; }
         for ( y_bin=0; y_bin<y_maxBin; y_bin++ ) { y_dens[y_bin] = 0; }
         for ( z_bin=0; z_bin<z_maxBin; z_bin++ ) { z_dens[z_bin] = 0; }

         for ( i=0; i<nAtoms; i++ ) 
         {
            x_bin = (int)( rx[i] * InvD ) ;
            y_bin = (int)( ry[i] * InvD ) ;
            z_bin = (int)( rz[i] * InvD ) ;
            x_dens[x_bin]++;
            y_dens[y_bin]++;
            z_dens[z_bin]++;
         }

         x_ref = 0;
         y_ref = 0;
         z_ref = 0;

         for ( x_bin=0; x_bin<x_maxBin; x_bin++ ) 
         { 
            if ( x_dens[x_bin] > x_ref )
            {
               x_ref = x_dens[x_bin]; 
               x_o = ((double)(x_bin)+0.5)*dR;
            }
         }
         for ( y_bin=0; y_bin<y_maxBin; y_bin++ ) 
         { 
            if ( y_dens[y_bin] > y_ref )
            {
               y_ref = y_dens[y_bin]; 
               y_o = ((double)(y_bin)+0.5)*dR;
            }
         }
         for ( z_bin=0; z_bin<z_maxBin; z_bin++ ) 
         { 
            if ( z_dens[z_bin] > z_ref )
            {
               z_ref = z_dens[z_bin]; 
               z_o = ((double)(z_bin)+0.5)*dR;
            }
         }


/*
 *       check the distance between droplet center and box center
 *       and shift the droplet to box center
 */
         x_o = 0.5*vec_box.x - x_o;
         y_o = 0.5*vec_box.y - y_o;
         z_o = 0.5*vec_box.z - z_o;
         for ( i=0; i<nAtoms; i++ )
         {
            rx[i] += x_o;
            ry[i] += y_o;
            rz[i] += z_o;
            if      ( rx[i] < 0.0 )       { rx[i] += vec_box.x; }
            else if ( rx[i] > vec_box.x ) { rx[i] -= vec_box.x; }
            if      ( ry[i] < 0.0 )       { ry[i] += vec_box.y; }
            else if ( ry[i] > vec_box.y ) { ry[i] -= vec_box.y; }
            if      ( rz[i] < 0.0 )       { rz[i] += vec_box.z; }
            else if ( rz[i] > vec_box.z ) { rz[i] -= vec_box.z; }
         }


/*
 *       find center of mass   
 */
         xcom = 0.0;
         ycom = 0.0;
         zcom = 0.0;
         msum = 0.0;
         for ( i=0; i<nAtoms; i++ ) 
         {
            mi = atom_param[iMol[i]][iAtom[i]].mass;
            xcom += mi*rx[i];
            ycom += mi*ry[i];
            zcom += mi*rz[i];
            msum += mi;
         }
         xcom /= msum;
         ycom /= msum;
         zcom /= msum;
/*
 *       center all atoms around COM and
 *       apply PBC to each molecule with respect to its central atom
 */
         for ( i=0; i<nAtoms; i++ ) 
         {
            rx[i] -= xcom;
            ry[i] -= ycom;
            rz[i] -= zcom;
         }

         for ( im=0; im<nMols; im++ )
         {
            iMin = mini[im][0];
            iMax = maxi[im][resNr[iMol[mini[im][0]]]-1];

            j = (iMin+iMax)/2;
            if      ( rx[j] < -0.5*vec_box.x ) { rx[j] += vec_box.x; }
            else if ( rx[j] >  0.5*vec_box.x ) { rx[j] -= vec_box.x; }
            if      ( ry[j] < -0.5*vec_box.y ) { ry[j] += vec_box.y; }
            else if ( ry[j] >  0.5*vec_box.y ) { ry[j] -= vec_box.y; }
            if      ( rz[j] < -0.5*vec_box.z ) { rz[j] += vec_box.z; }
            else if ( rz[j] >  0.5*vec_box.z ) { rz[j] -= vec_box.z; }

            for (i=iMin; i<iMax; i++) 
            {
               if ( i==j ) { continue; }
               if      ( rx[i] - rx[j] < -0.5*vec_box.x ) { rx[i] += vec_box.x; }
               else if ( rx[i] - rx[j] >  0.5*vec_box.x ) { rx[i] -= vec_box.x; }
               if      ( ry[i] - ry[j] < -0.5*vec_box.y ) { ry[i] += vec_box.y; }
               else if ( ry[i] - ry[j] >  0.5*vec_box.y ) { ry[i] -= vec_box.y; }
               if      ( rz[i] - rz[j] < -0.5*vec_box.z ) { rz[i] += vec_box.z; }
               else if ( rz[i] - rz[j] >  0.5*vec_box.z ) { rz[i] -= vec_box.z; }
            }
         }

         for ( i=0; i<nAtoms; i++ ) 
         {
            r_ci[i] = sqrt(rx[i]*rx[i] + ry[i]*ry[i] + rz[i]*rz[i]);
         }

/*
 *       initialize arrays: dens and pU
 */
         for ( bin=0; bin<maxBin; bin++ ) 
         {
            dens[bin] = 0.0;
            for ( mol=0; mol<molTypes; mol++ ) 
            {
               for ( res=0; res<resNr[mol]; res++ )
               {
                  resDens[bin][mol][res] = 0.0;
               }
            }

            pN_U_LJ[bin]  = 0.0;
            pN_U_LJc[bin] = 0.0;
            pN_U_QQ[bin]  = 0.0;
            pN_U_14[bin]  = 0.0;
            pN_U_BD1[bin] = 0.0;
            pN_U_BD2[bin] = 0.0;
            pN_K[bin] = 0.0;
            pN[bin]  = 0.0;
            pNc[bin] = 0.0;

            pT_U_LJ[bin]  = 0.0;
            pT_U_LJc[bin] = 0.0;
            pT_U_QQ[bin]  = 0.0;
            pT_U_14[bin]  = 0.0;
            pT_U_BD1[bin] = 0.0;
            pT_U_BD2[bin] = 0.0;
            pT_K[bin] = 0.0;
            pT[bin]  = 0.0;
            pTc[bin] = 0.0;
         }

/*
 *       generate density profile   
 */
         im = -1;
         for ( mol=0; mol<molTypes; mol++ )
         {
            for ( km=0; km<molNr[mol]; km++ )
            {
               im++;
               for ( res=0; res<resNr[mol]; res++ )
               {
/*
                  cx[im][res] = 0.0;
                  cy[im][res] = 0.0;
                  cz[im][res] = 0.0;
                  msum = 0.0;
*/
                  iMin = mini[im][res];
                  iMax = maxi[im][res];
                  for (i=iMin; i<iMax; i++) 
                  {
                     vec_ri.x = rx[i];
                     vec_ri.y = ry[i];
                     vec_ri.z = rz[i];
                     rci = r_ci[i];
                     bin = (int)( rci * InvD );
   
                     vec_vi.x = vx[i];
                     vec_vi.y = vy[i];
                     vec_vi.z = vz[i];
                     vci2 = dist_2( vec_vi );
   
                     mi = atom_param[iMol[i]][iAtom[i]].mass;
                     kE_N = mi*pow(inner_product(vec_vi,vec_ri),2.0)/(rci*rci);
                     kE_T = 0.5*(mi*vci2 - kE_N);
   
                     if ( bin < maxBin )
                     {
                        incr_dens = mi * inv_dV[bin]; /* unit: g mol^-1 nm^-3 */
                        incr_dens /= nA*100.0; /* unit: g cm^-3 */
                        dens[bin] += incr_dens;
                        resDens[bin][iMol[iMin]][res] += incr_dens;

                        pN_K[bin] += kE_N * inv_dV[bin];

                        pT_K[bin] += kE_T * inv_dV[bin];
                     }
                  }
               }
            }
         }

/*
 *       loop over all molecules to calculate pU
 *       from intra-mol and inter-mol contributions
 */

/*
 *       bonded interaction
 */
         for ( im=0; im<nMols; im++ )
         {

            mol = iMol[mini[im][0]];

/*
 *          bonds
 */
            for ( iBond=0; iBond<nBonds[mol]; iBond++ )
            {
               if ( bond_funct[mol][iBond] == 1 )
               {
                  i = bond_1[mol][iBond].atom_i + mini[im][0];
                  j = bond_1[mol][iBond].atom_j + mini[im][0];

                  vec_ri.x = rx[i];
                  vec_ri.y = ry[i];
                  vec_ri.z = rz[i];
                  rci = r_ci[i];

                  vec_rj.x = rx[j];
                  vec_rj.y = ry[j];
                  vec_rj.z = rz[j];
                  rcj = r_ci[j];

                  b0 = bond_1[mol][iBond].b0;
                  kb = bond_1[mol][iBond].kb;
                  bond_force = compute_bond_1 (vec_ri.x, vec_ri.y, vec_ri.z,
                                               vec_rj.x, vec_rj.y, vec_rj.z,
                                               b0, kb);

                  compute_pres_BD_1 ( maxBin, dR, InvD, inv_dV, thrshd,
                                      vec_ri, rci, vec_rj, rcj, bond_force.fi );
               }
               else if ( bond_funct[mol][iBond] == 3 )
               {
                  i = bond_3[mol][iBond].atom_i + mini[im][0];
                  j = bond_3[mol][iBond].atom_j + mini[im][0];

                  vec_ri.x = rx[i];
                  vec_ri.y = ry[i];
                  vec_ri.z = rz[i];
                  rci = r_ci[i];

                  vec_rj.x = rx[j];
                  vec_rj.y = ry[j];
                  vec_rj.z = rz[j];
                  rcj = r_ci[j];

                  b0 = bond_3[mol][iBond].b0;
                  D = bond_3[mol][iBond].D;
                  beta = bond_3[mol][iBond].beta;
                  bond_force = compute_bond_3 (vec_ri.x, vec_ri.y, vec_ri.z,
                                               vec_rj.x, vec_rj.y, vec_rj.z,
                                               b0, D, beta);

                  compute_pres_BD_1 ( maxBin, dR, InvD, inv_dV, thrshd,
                                      vec_ri, rci, vec_rj, rcj, bond_force.fi );
               }
            }

/*
 *          pairs
 */
            for ( iPair=0; iPair<nPairs[mol]; iPair++ )
            {
               if ( pair_funct[mol][iPair] == 1 )
               {
                  i = pair_1[mol][iPair].atom_i + mini[im][0];
                  j = pair_1[mol][iPair].atom_j + mini[im][0];

                  qi    = atom_param[iMol[i]][iAtom[i]].charge;
                  iType = atom_param[iMol[i]][iAtom[i]].atomtype;
                  qj    = atom_param[iMol[j]][iAtom[j]].charge;
                  jType = atom_param[iMol[j]][iAtom[j]].atomtype;
                  c6   =   LJ_C6[iType][jType];
                  c12  =  LJ_C12[iType][jType];
/*
 *                Note i and j are interchanged in compute_LJ_QQ
 *                so that the force on i is computed
 */
                  LJ_QQ = compute_LJ_QQ (rx[j],ry[j],rz[j],
                                         rx[i],ry[i],rz[i],
                                         fudgeLJ,fudgeQQ,
                                         c6,c12,rCut,fQQ,qi,qj);
               }
               else if ( pair_funct[mol][iPair] == 2 )
               {
                  i = pair_2[mol][iPair].atom_i + mini[im][0];
                  j = pair_2[mol][iPair].atom_j + mini[im][0];

                  pair_fudgeQQ = pair_2[mol][iPair].fudgeQQ;
                  qi   =  pair_2[mol][iPair].q_i;
                  qj   =  pair_2[mol][iPair].q_j;
/*          
 *                V = sigma, W = epsilon for comb_rule = 2 or 3 ( 1 not supported )
 */
                  c6   =  pair_2[mol][iPair].W * pow(pair_2[mol][iPair].V,6.0) * 4.0;
                  c12  =  pair_2[mol][iPair].W * pow(pair_2[mol][iPair].V,12.0) * 4.0;
/*
 *                Note i and j are interchanged in compute_LJ_QQ
 *                so that the force on i is computed
 */
                  LJ_QQ = compute_LJ_QQ (rx[j],ry[j],rz[j],
                                         rx[i],ry[i],rz[i],
                                         1.0,pair_fudgeQQ,
                                         c6,c12,rCut,fQQ,qi,qj);
               }

               compute_pres_14_1 ( maxBin, dR, InvD, inv_dV, thrshd,
                                   vec_ri, rci, vec_rj, rcj, 
                                   add_vector(LJ_QQ.LJ,LJ_QQ.QQ) );
            }

/*
 *          angles
 */
            for ( iAngle=0; iAngle<nAngles[mol]; iAngle++ )
            {
               if ( angle_funct[mol][iAngle] == 1 )
               {
                  i = angle_1[mol][iAngle].atom_i + mini[im][0];
                  j = angle_1[mol][iAngle].atom_j + mini[im][0];
                  k = angle_1[mol][iAngle].atom_k + mini[im][0];

                  a0 = angle_1[mol][iAngle].a0;
                  ka = angle_1[mol][iAngle].ka;
                  angle_force = compute_angle_1 ( rx[i], ry[i], rz[i],
                                                  rx[j], ry[j], rz[j],
                                                  rx[k], ry[k], rz[k],
                                                  a0, ka );
                   
                  vec_rj.x = rx[j];
                  vec_rj.y = ry[j];
                  vec_rj.z = rz[j];
                  rcj = r_ci[j];

                  for ( iia=0; iia<2; iia++ )
                  {
                     if ( iia == 0 )
                     {
                        vec_ri.x = rx[i];
                        vec_ri.y = ry[i];
                        vec_ri.z = rz[i];
                        rci = r_ci[i];
                        compute_pres_BD_2 ( maxBin, dR, InvD, inv_dV, thrshd,
                                            vec_ri, rci, vec_rj, rcj, angle_force.fi );
                     }
                     else if ( iia == 1 )
                     {
                        vec_ri.x = rx[k];
                        vec_ri.y = ry[k];
                        vec_ri.z = rz[k];
                        rci = r_ci[k];
                        compute_pres_BD_2 ( maxBin, dR, InvD, inv_dV, thrshd,
                                            vec_ri, rci, vec_rj, rcj, angle_force.fk );
                     }
                  }
               }
            }

/*
 *          dihedrals
 */
            for ( iDihedral=0; iDihedral<nDihedrals[mol]; iDihedral++ )
            {
               if ( dihedral_funct[mol][iDihedral] == 3 )
               {
                     i = dihedral_3[mol][iDihedral].atom_i + mini[im][0];
                     j = dihedral_3[mol][iDihedral].atom_j + mini[im][0];
                     k = dihedral_3[mol][iDihedral].atom_k + mini[im][0];
                     l = dihedral_3[mol][iDihedral].atom_l + mini[im][0];

                     c0 = dihedral_3[mol][iDihedral].c0;
                     c1 = dihedral_3[mol][iDihedral].c1;
                     c2 = dihedral_3[mol][iDihedral].c2;
                     c3 = dihedral_3[mol][iDihedral].c3;
                     c4 = dihedral_3[mol][iDihedral].c4;
                     c5 = dihedral_3[mol][iDihedral].c5;

                     dihedral_force = compute_dihedral_3 ( rx[i], ry[i], rz[i],
                                                           rx[j], ry[j], rz[j],
                                                           rx[k], ry[k], rz[k],
                                                           rx[l], ry[l], rz[l],
                                                           c0, c1, c2, c3, c4, c5 );
               }
               else if ( dihedral_funct[mol][iDihedral]==1 ||
                         dihedral_funct[mol][iDihedral]==4 ||
                         dihedral_funct[mol][iDihedral]==9 )
               {
                     i = dihedral_9[mol][iDihedral].atom_i + mini[im][0];
                     j = dihedral_9[mol][iDihedral].atom_j + mini[im][0];
                     k = dihedral_9[mol][iDihedral].atom_k + mini[im][0];
                     l = dihedral_9[mol][iDihedral].atom_l + mini[im][0];

                     phi0 = dihedral_9[mol][iDihedral].phi0;
                     kphi = dihedral_9[mol][iDihedral].kphi;
                     n = dihedral_9[mol][iDihedral].n;

                     dihedral_force = compute_dihedral_9 ( rx[i], ry[i], rz[i],
                                                           rx[j], ry[j], rz[j],
                                                           rx[k], ry[k], rz[k],
                                                           rx[l], ry[l], rz[l],
                                                           phi0, kphi, n );
               }

               vec_rj.x = rx[j];
               vec_rj.y = ry[j];
               vec_rj.z = rz[j];
               rcj = r_ci[j];

               for ( iia=0; iia<3; iia++ )
               {
                  if ( iia == 0 )
                  {
                     vec_ri.x = rx[i];
                     vec_ri.y = ry[i];
                     vec_ri.z = rz[i];
                     rci = r_ci[i];
                     compute_pres_BD_2 ( maxBin, dR, InvD, inv_dV, thrshd,
                                         vec_ri, rci, vec_rj, rcj, dihedral_force.fi );
                  }
                  else if ( iia == 1 )
                  {
                     vec_ri.x = rx[k];
                     vec_ri.y = ry[k];
                     vec_ri.z = rz[k];
                     rci = r_ci[k];
                     compute_pres_BD_2 ( maxBin, dR, InvD, inv_dV, thrshd,
                                         vec_ri, rci, vec_rj, rcj, dihedral_force.fk );
                  }
                  else if ( iia == 2 )
                  {
/*
 *                   Note letter l is not number 1
 */
                     vec_ri.x = rx[l];
                     vec_ri.y = ry[l];
                     vec_ri.z = rz[l];
                     rci = r_ci[l];
                     compute_pres_BD_2 ( maxBin, dR, InvD, inv_dV, thrshd,
                                         vec_ri, rci, vec_rj, rcj, dihedral_force.fl );
                  }
               }
            }

         }

/*
 *       intramolecular nonbonded interaction
 */
         for ( im=0; im<nMols; im++ )
         {
            iMin = mini[im][0];
            iMax = maxi[im][resNr[iMol[mini[im][0]]]-1];

            for (i=iMin; i<iMax; i++) 
            {
            for (j=i+1; j<iMax; j++)
            {

               if ( exclude[iMol[mini[im][0]]][iAtom[i]][iAtom[j]]==1 )
               {
                  continue;
               }

               vec_ri.x = rx[i];
               vec_ri.y = ry[i];
               vec_ri.z = rz[i];
               rci = r_ci[i];

               vec_rj.x = rx[j];
               vec_rj.y = ry[j];
               vec_rj.z = rz[j];
               rcj = r_ci[j];

               qi    = atom_param[iMol[i]][iAtom[i]].charge;
               iType = atom_param[iMol[i]][iAtom[i]].atomtype;
               qj    = atom_param[iMol[j]][iAtom[j]].charge;
               jType = atom_param[iMol[j]][iAtom[j]].atomtype;
               c6   =   LJ_C6[iType][jType];
               c12  =  LJ_C12[iType][jType];
/*
 *             Note i and j are interchanged in compute_LJ_QQ
 *             so that the force on i is computed
 */
               LJ_QQ = compute_LJ_QQ (rx[j],ry[j],rz[j],
                                      rx[i],ry[i],rz[i],
                                      1.0,1.0,
                                      c6,c12,rCut,fQQ,qi,qj);

               compute_pres_NB_1 ( maxBin, dR, InvD, inv_dV, thrshd,
                                   vec_ri, rci, vec_rj, rcj, LJ_QQ );

            }
            }

         }

/*
 *       intermolecular nonbonded interaction
 */
         for ( im=0; im<nMols; im++ )
         {
            iMin = mini[im][0];
            iMax = maxi[im][resNr[iMol[mini[im][0]]]-1];

         for ( jm=im+1; jm<nMols; jm++ )
         {
            jMin = mini[jm][0];
            jMax = maxi[jm][resNr[iMol[mini[jm][0]]]-1];

            for (i=iMin; i<iMax; i++) 
            {
            for (j=jMin; j<jMax; j++) 
            {

               vec_ri.x = rx[i];
               vec_ri.y = ry[i];
               vec_ri.z = rz[i];
               rci = r_ci[i];

               vec_rj.x = rx[j];
               vec_rj.y = ry[j];
               vec_rj.z = rz[j];
               rcj = r_ci[j];

               qi    = atom_param[iMol[i]][iAtom[i]].charge;
               iType = atom_param[iMol[i]][iAtom[i]].atomtype;
               qj    = atom_param[iMol[j]][iAtom[j]].charge;
               jType = atom_param[iMol[j]][iAtom[j]].atomtype;
               c6   =   LJ_C6[iType][jType];
               c12  =  LJ_C12[iType][jType];
/*
 *             Note i and j are interchanged in compute_LJ_QQ
 *             so that the force on i is computed
 */
               LJ_QQ = compute_LJ_QQ (rx[j],ry[j],rz[j],
                                      rx[i],ry[i],rz[i],
                                      1.0,1.0,
                                      c6,c12,rCut,fQQ,qi,qj);

               compute_pres_NB_1 ( maxBin, dR, InvD, inv_dV, thrshd,
                                   vec_ri, rci, vec_rj, rcj, LJ_QQ );

            }
            }

         }
         }

/*
 *       calculate work of formation  
 */
         work_N_full = 0.0;
         work_N_cut  = 0.0;
         work_T_full = 0.0;
         work_T_cut  = 0.0;
         for ( bin=0; bin<maxBin; bin++ ) 
         {
            r = ( (double)(bin) + 0.5 ) * dR;

            pN[bin]  = pN_K[bin]+pN_U_LJ[bin]+pN_U_QQ[bin]+pN_U_14[bin]+pN_U_BD1[bin]+pN_U_BD2[bin];
            pNc[bin] = pN_K[bin]+pN_U_LJc[bin]+pN_U_QQ[bin]+pN_U_14[bin]+pN_U_BD1[bin]+pN_U_BD2[bin];

            work_N_full += pN[bin]*r*r * dR;
            work_N_cut += pNc[bin]*r*r * dR;

            pT[bin]  = pT_K[bin]+pT_U_LJ[bin]+pT_U_QQ[bin]+pT_U_14[bin]+pT_U_BD1[bin]+pT_U_BD2[bin];
            pTc[bin] = pT_K[bin]+pT_U_LJc[bin]+pT_U_QQ[bin]+pT_U_14[bin]+pT_U_BD1[bin]+pT_U_BD2[bin];

            work_T_full += pT[bin]*r*r * dR;
            work_T_cut += pTc[bin]*r*r * dR;
         }
         work_N_full *=  pi*2.0;
         work_N_cut  *=  pi*2.0;
         work_T_full *= -pi*4.0;
         work_T_cut  *= -pi*4.0;

         printf ( "%s\n", "#Frame  W_N_full  W_N_cut  W_T_full  W_T_cut" );
         printf ( "%d  %.6f  %.6f  %.6f  %.6f\n", 
                  frame, work_N_full, work_N_cut, work_T_full, work_T_cut );


/*
 *       write results
 */ 
         printf ( "%s\n", "#Radius  pN  pNc  pN_K  pN_U_LJ  pN_U_LJc  pN_U_QQ  pN_U_14  pN_U_BD1  pN_U_BD2" );
         for ( bin=0; bin<maxBin; bin++ ) 
         {
            r = ( (double)(bin) + 0.5 ) * dR;
            printf ( "%.3f  %.6f %.6f %.6f  %.6f %.6f %.6f  %.6f %.6f %.6f\n", 
                     r, pN[bin], pNc[bin], pN_K[bin],
                     pN_U_LJ[bin], pN_U_LJc[bin], pN_U_QQ[bin],
                     pN_U_14[bin], pN_U_BD1[bin], pN_U_BD2[bin] );
         }
       
         printf ( "%s\n", "#Radius  pT  pTc  pT_K  pT_U_LJ  pT_U_LJc  pT_U_QQ  pT_U_14  pT_U_BD1  pT_U_BD2" );
         for ( bin=0; bin<maxBin; bin++ ) 
         {
            r = ( (double)(bin) + 0.5 ) * dR;
            printf ( "%.3f  %.6f %.6f %.6f  %.6f %.6f %.6f  %.6f %.6f %.6f\n", 
                     r, pT[bin], pTc[bin], pT_K[bin],
                     pT_U_LJ[bin], pT_U_LJc[bin], pT_U_QQ[bin],
                     pT_U_14[bin], pT_U_BD1[bin], pT_U_BD2[bin] );
         }
/*     
 *       write molecular densities
 */ 
         printf ( "%s\n", "#Radius  Total_Density  Residue_Densities (g/cm^3)" );
         for ( bin=0; bin<maxBin; bin++ ) 
         {
            r = ( (double)(bin) + 0.5 ) * dR;
            printf ( "%.3f  %.6f", r, dens[bin] );
            for ( mol=0; mol<molTypes; mol++ )
            {
               for ( res=0; res<resNr[mol]; res++ )
               {
                  printf ( "  %.6f", resDens[bin][mol][res] );
               }
            }
            printf ( "\n" );
         }

/*
 *    end of loop (frame)   
 */
      }
/*
 *    close trajectory file
 */
      xdrfile_close (trr);


/*****************************************************************************
 *    End of function main                                                   *
 *****************************************************************************/

      end_t = time(NULL);
      delta_time = (int)(difftime(end_t,start_t));

      printf ( ">>> Calculation ended normally\n" );
      printf ( ">>> %d seconds were used\n", delta_time );

   return 0;
   }

