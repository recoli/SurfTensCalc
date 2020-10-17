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
 *              Version  1.6-test, 2012-01-07                       *
 ********************************************************************/

   #include <math.h>
   #include <stdio.h>
   #include <stdlib.h>
   #include <xdrfile_xtc.h>

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
         printf("%s\n", "Incorrect number of arguments (should be 5)!");
         printf("%s\n", "1st argument: nFrame");
         printf("%s\n", "2nd argument: nStart");
         printf("%s\n", "3rd argument: temperature (in Kelvin)");
         printf("%s\n", "4th argument: cutoff radius (in nm)");
         printf("%s\n", "5th argument: approx. box size (in nm)");
         exit(1);
      }
      printf ( "%s\n", "#Frame  W_Full  W_Cut" );

/********************************************************************
 *    Define variables                                              *
 ********************************************************************/

/*
 *    number of frames in trajectory file 
 */
      int   frame, nFrames, nStart;
      nFrames = atoi(argv[1]);
      nStart = atoi(argv[2]);

      if ( nFrames<=nStart ) 
      {
         printf( "Error: nFrames is no larger than nStart!\n" ) ;
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

/*
      int comb_rule;
      double fudgeLJ, fudgeQQ;
*/

      int **atomtype, **residue;

      int nTypes;
      double **LJ_C6, **LJ_C12;
      int iType, jType;
      double c6, c12;

/*
 *    work of formation
 */
      double workFull, workCut;
/*
 *    coordinates of atoms and COMs, COM = center of mass 
 */
      double *rx, *ry, *rz;
      double **cx, **cy, **cz;
/*
 *    thickness of sperical layer, 0.3 Angstrom = 0.1 Sigma_Oxygen
 */  
      const double dR = 0.03 ;
      int maxBin ;
      maxBin = (int)( boxSize * 0.4 / dR ) ;

      int    bin ;
      double dV ;
      double *averDens, *dens ;
      double *averPULJ, *averPUQQ, *averPULJc;
      double *pULJ, *pUQQ, *pULJc;
      double pK;
/*
 *    densities for different molecules 
 */
      double **molDens ;
      double ***resDens ;
/*
 *    increment, temporary variable 
 */
      double incr, incrLJ, incrQQ, incrLJc;
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
      int    i, j, inters, im, jm, iMin, iMax, jMin, jMax ;
      double mi, qi, qj ;
      double r, r2 ;
      double chord2 ;

      double **r_ci2, **r_ci;

      double cxi, cyi, czi, rci2, rci ;
      double cxj, cyj, czj, rcj2, rcj ;
      double cxij, cyij, czij, cij2, cij ; 
      double rxij, ryij, rzij, rij2, rij ;
      double rij6, rij12, fij ;
      double ffijLJ, ffijQQ, ffijLJc;
      double fxijLJ, fxijQQ, fxijLJc;
      double fyijLJ, fyijQQ, fyijLJc;
      double fzijLJ, fzijQQ, fzijLJc;
/*
 *    variables for interception calculation   
 */
      double cosine, cos_i, h2 ;
      double rMin, rMax ;
      int    bMin, bMax ;
/*
 *    variables for xdrfile_xtc  
 */
      int     step,read_return,gmxNAtoms;
      float   time,prec;
      matrix  box;
      rvec    *x;
      XDRFILE *xtc;



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
      char    line[256], tmp[50], filenameXTC[50];
      double  **mass, **charge, **sigma, **epsilon;
      int     nAtoms, nMols;
      double  *molMass;
      double  **resMass;
/*
 *    open param.txt 
 */
      file_par = fopen("param.txt","r") ;
      if ( NULL==file_par ) {
         printf( "Cannot open file: param.txt !\n" ) ;
         exit(1);
      }
      mol = 0;
/*
 *    read comb_rule and fudgeLJ, fudgeQQ (first two lines)
 */
      fgets( line, sizeof( line ), file_par );
/*
      sscanf( line, "%s%d", tmp, &comb_rule );
*/
      fgets( line, sizeof( line ), file_par );
/*
      sscanf( line, "%s%f%f", tmp, &fudgeLJ, &fudgeQQ );
*/
/*
 *    read name of xtc file
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%s", tmp, filenameXTC );
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
      charge  = malloc(molTypes * sizeof(double *));
      mass    = malloc(molTypes * sizeof(double *));
      sigma   = malloc(molTypes * sizeof(double *));
      epsilon = malloc(molTypes * sizeof(double *));
      atomtype = malloc(molTypes * sizeof(int *));
      residue = malloc(molTypes * sizeof(int *));
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &atom_in_mol[mol] );
         charge[mol]  = malloc(atom_in_mol[mol] * sizeof(double));
         mass[mol]    = malloc(atom_in_mol[mol] * sizeof(double));
         sigma[mol]   = malloc(atom_in_mol[mol] * sizeof(double));
         epsilon[mol] = malloc(atom_in_mol[mol] * sizeof(double));
         atomtype[mol] = malloc(atom_in_mol[mol] * sizeof(int));
         residue[mol] = malloc(atom_in_mol[mol] * sizeof(int));
         for ( atom=0; atom<atom_in_mol[mol]; atom++ ) 
         {
            fgets( line, sizeof( line ), file_par );
            sscanf( line, "%lf%lf%lf%lf%d%d",
                    &charge[mol][atom], &mass[mol][atom],
                    &sigma[mol][atom],  &epsilon[mol][atom],
                    &atomtype[mol][atom], &residue[mol][atom] );
            res = residue[mol][atom];
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
            molMass[mol] += mass[mol][atom];
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
               resMass[mol][res] += mass[mol][atom];
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
        printf("Incorrect number of atoms in the system!\n");
        printf("Program stops at assigning atom index and mol index!\n");
        exit(1);
      }

      if ( nMols != im ) 
      {
        printf("Incorrect number of molecules in the system!\n");
        printf("Program stops at assigning atom index and mol index!\n");
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
        printf("Incorrect number of atoms in the system!\n");
        printf("Program stops at assigning atom index and mol index!\n");
        exit(1);
      }

      if ( nMols != im ) 
      {
        printf("Incorrect number of molecules in the system!\n");
        printf("Program stops at assigning atom index and mol index!\n");
        exit(1);
      }

/********************************************************************
 *    Allocating arrays                                             *
 ********************************************************************/

      rx = malloc( sizeof(double) * nAtoms );
      ry = malloc( sizeof(double) * nAtoms );
      rz = malloc( sizeof(double) * nAtoms );
      if ( NULL==rx || NULL==ry || NULL==rz ) 
      {
         printf ("Unable to allocate arrays: rx, ry, rz!\n");
         exit(1);
      }

      cx = malloc( sizeof(double *) * nMols );
      cy = malloc( sizeof(double *) * nMols );
      cz = malloc( sizeof(double *) * nMols );
      r_ci2 = malloc( sizeof(double *) * nMols );
      r_ci  = malloc( sizeof(double *) * nMols );
      im = 0;
      for ( mol=0; mol<molTypes; mol++ )
      {
         for ( km=0; km<molNr[mol]; km++ )
         {
            cx[im]  = malloc(resNr[mol] * sizeof(double));
            cy[im]  = malloc(resNr[mol] * sizeof(double));
            cz[im]  = malloc(resNr[mol] * sizeof(double));
            r_ci2[im]  = malloc(resNr[mol] * sizeof(double));
            r_ci[im]   = malloc(resNr[mol] * sizeof(double));
            im++;
         }
      }
      if ( NULL==cx || NULL==cy || NULL==cz || NULL==r_ci2 || NULL==r_ci ) 
      {
         printf ("Unable to allocate arrays: cx, cy, cz, r_ci2, r_ci!\n");
         exit(1);
      }

      averDens  = malloc( sizeof( double ) * maxBin );
      averPULJ  = malloc( sizeof( double ) * maxBin );
      averPUQQ  = malloc( sizeof( double ) * maxBin );
      averPULJc = malloc( sizeof( double ) * maxBin );
      if ( NULL==averDens || NULL==averPULJ || NULL==averPUQQ || NULL==averPULJc ) 
      {
         printf ("Unable to allocate arrays: averDens, averPULJ, averPUQQ, averPULJc!\n");
         exit(1);
      }

      dens  = malloc( sizeof( double ) * maxBin );
      pULJ  = malloc( sizeof( double ) * maxBin );
      pUQQ  = malloc( sizeof( double ) * maxBin );
      pULJc = malloc( sizeof( double ) * maxBin );
      if ( NULL==dens || NULL==pULJ || NULL==pUQQ || NULL==pULJc ) 
      {
         printf ("Unable to allocate arrays: dens, pULJ, pUQQ, pULJc!\n");
         exit(1);
      }

      molDens = malloc(maxBin * sizeof(double *));
      resDens = malloc(maxBin * sizeof(double **));
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         molDens[bin] = malloc(molTypes * sizeof(double));
         resDens[bin] = malloc(molTypes * sizeof(double *));
         for ( mol=0; mol<molTypes; mol++ )
         {
            resDens[bin][mol] = malloc(resNr[mol] * sizeof(double));
         }
      }

      x = calloc(nAtoms,sizeof(x[0]));
      if ( NULL==x ) 
      {
         printf ("Unable to allocate array: x!\n");
         exit(1);
      }

      totalPairs = nMols*(nMols-1)/2 ;
      molI = malloc( sizeof( int ) * totalPairs );
      molJ = malloc( sizeof( int ) * totalPairs );
      if ( NULL==molI || NULL==molJ ) 
      {
         printf ("Unable to allocate arrays: molI, molJ!\n");
         exit(1);
      }


/********************************************************************
 *    Calculating pressure profile and work of formation            *
 ********************************************************************/

/*
 *    initialize arrays: averDens, averPU, molDens
 */
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         averDens[bin]  = 0.0;
         averPULJ[bin]  = 0.0;
         averPUQQ[bin]  = 0.0;
         averPULJc[bin] = 0.0;
         for ( mol=0; mol<molTypes; mol++ ) 
         {
            molDens[bin][mol] = 0.0;
            for ( res=0; res<resNr[mol]; res++ )
            {
               resDens[bin][mol][res] = 0.0;
            }
         }
      }
/*
 *    read xtc files for coordinates   
 */
      xtc = xdrfile_open (filenameXTC,"r");
      read_xtc_natoms (filenameXTC, &gmxNAtoms);
      /*printf("Total Number of Frames %10d\n", nFrames);*/
      /*printf("Starting from Frame    %10d\n", nStart);*/
/*
 *    start of loop (frame)   
 */
      for (frame=0; frame<nFrames; frame++) 
      {
/*
 *       read xtc file   
 */
         read_return = read_xtc (xtc,gmxNAtoms,&step,&time,box,x,&prec);
         if (read_return!=0) {
            printf ("Failed to read xtc file!\n");
            exit(1);
         }
         if (gmxNAtoms!=nAtoms) {
            printf ("Incorrect number of Atoms!\n");
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
/*
 *          mi = mass for atom i
 */
            mi = mass[iMol[i]][iAtom[i]];
            xcom += mi*rx[i];
            ycom += mi*ry[i];
            zcom += mi*rz[i];
            msum += mi;
         }
         xcom /= msum;
         ycom /= msum;
         zcom /= msum;
/*
 *       center all atoms around COM   
 */
         for ( i=0; i<nAtoms; i++ ) 
         {
            rx[i] -= xcom;
            ry[i] -= ycom;
            rz[i] -= zcom;
         }
/*
 *       initialize arrays: dens and pU
 */
         for ( bin=0; bin<maxBin; bin++ ) 
         {
            dens[bin] = 0.0;
            pULJ[bin] = 0.0;
            pUQQ[bin] = 0.0;
            pULJc[bin]= 0.0;
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
                  cx[im][res] = 0.0;
                  cy[im][res] = 0.0;
                  cz[im][res] = 0.0;
                  msum = 0.0;
                  iMin = mini[im][res];
                  iMax = maxi[im][res];
                  for (i=iMin; i<iMax; i++) 
                  {
                     mi = mass[iMol[i]][iAtom[i]];
                     cx[im][res] += mi*rx[i];
                     cy[im][res] += mi*ry[i];
                     cz[im][res] += mi*rz[i];
                     msum += mi;
                  }
                  if ( fabs(msum-resMass[iMol[iMin]][res])>0.01 ) 
                  {
                     printf ("Incorrect molar weight of molecule %d residue %d!\n", im, res+1);
                     printf ("msum = %f, resMass[%d][%d]=%f!\n", 
                             msum, iMol[iMin], res, resMass[iMol[iMin]][res]);
                     exit(1);
                  }
                  cx[im][res] /= msum;
                  cy[im][res] /= msum;
                  cz[im][res] /= msum;

                  r_ci2[im][res] = cx[im][res]*cx[im][res] 
                                 + cy[im][res]*cy[im][res] 
                                 + cz[im][res]*cz[im][res];
                  r_ci[im][res] = sqrt( r_ci2[im][res] );

                  bin = (int)( r_ci[im][res]/dR );

                  if ( bin < maxBin ) 
                  {
/*            
 *                   volume of the spherical layer 
 */           
                     dV = pi*4.0/3.0*( pow((double)(bin+1)*dR,3.0) - pow((double)bin*dR,3.0) );
/*            
 *                   add increment to dens[bin] and averDens[bin]
 */           
                     incr = 1.0/dV;
                     dens[bin] += incr;
                     averDens[bin] += incr;
                     resDens[bin][iMol[iMin]][res] += incr;
                  }
               }
            }
         }
/*
 *       find all pairs   
 */
         pair = 0;
         for ( im=0; im<nMols-1; im++ ) 
         {
            for ( jm=im+1; jm<nMols; jm++ ) 
            {
               molI[pair] = im;
               molJ[pair] = jm;
               pair++;
            }
         }
         if ( totalPairs != pair ) 
         {
            printf ( "Incorrect number of pairs of molecules!\n" );
            exit(1);
         }
/*
 *       loop over all pairs to calculate pU
 */
         for ( pair=0; pair<totalPairs; pair++ ) 
         {
            im = molI[pair];
            jm = molJ[pair];

            for ( ires=0; ires<resNr[iMol[mini[im][0]]]; ires++ )
            {
            for ( jres=0; jres<resNr[iMol[mini[jm][0]]]; jres++ )
            {
/*
 *             center of mass distances   
 */          
               rci2 = r_ci2[im][ires];
               rci = r_ci[im][ires];
             
               rcj2 = r_ci2[jm][jres];
               rcj = r_ci[jm][jres];
             
               cxij = cx[jm][jres] - cx[im][ires];
               cyij = cy[jm][jres] - cy[im][ires];
               czij = cz[jm][jres] - cz[im][ires];
               cij2 = cxij*cxij + cyij*cyij + czij*czij;
               cij  = sqrt( cij2 );
/*           
 *             calculate forces between im and jm
 */          
               fxijLJ = 0.0;
               fyijLJ = 0.0;
               fzijLJ = 0.0;
             
               fxijQQ = 0.0;
               fyijQQ = 0.0;
               fzijQQ = 0.0;
             
               fxijLJc = 0.0;
               fyijLJc = 0.0;
               fzijLJc = 0.0;
/*           
 *             determine the upper and lower boundary of i and j   
 */          
               iMin = mini[im][ires];
               iMax = maxi[im][ires];
               jMin = mini[jm][jres];
               jMax = maxi[jm][jres];
/*           
 *             start of loops (i and j)   
 */          
               for (i=iMin; i<iMax; i++) 
               {
               for (j=jMin; j<jMax; j++) 
               {
/*           
 *                interatomic distances   
 */          
                  rxij = rx[j] - rx[i];
                  ryij = ry[j] - ry[i];
                  rzij = rz[j] - rz[i];
                  rij2 = rxij*rxij + ryij*ryij + rzij*rzij;
                  rij  = sqrt( rij2 );
/*           
 *                parameters for atom i and j
 */          
                  mol  = iMol[i];
                  atom = iAtom[i];
                  qi   =  charge[mol][atom];
                  iType = atomtype[mol][atom];
             
                  mol  = iMol[j];
                  atom = iAtom[j];
                  qj   =  charge[mol][atom];
                  jType = atomtype[mol][atom];
             
                  c6   =   LJ_C6[iType][jType];
                  c12  =  LJ_C12[iType][jType];
/*           
 *                VDW forces, LJ potential ( sigma, epsilon )
 *           
 *                U = C12 / R^12 - C6 / R^6
 *                F = 12 * C12 / R^13 - 6 * C6 / R^7
 *                vecF = ( 12 * C12 / R^14 - 6 * C6 / R^8 ) * vecR
 *                     = ( 2 * C12 / R^12 - C6 / R^6 ) * 6 / R^2 * vecR
 *           
 *                U = 4 * eps * ( sig^12/R^12 - sig^6/R^6 )
 *                F = 4 * eps * ( 12*sig^12/R^13 - 6*sig^6/R^7 )
 *                vecF = 4 * eps * ( 12*sig^12/R^14 - 6*sig^6/R^8 ) * vecR
 *                     = 24 * eps * ( 2*(sig/R)^12 - (sig/R)^6 ) / R^2 * vecR
 */          
                  if ( c6!=0.0 && c12!=0.0 ) 
                  {
                     rij6 = rij2 * rij2 * rij2 ;
                     rij12 = rij6 * rij6;
                     fij = ( c12*2.0/rij12 - c6/rij6 ) * 6.0/rij2;
                     fxijLJ = fxijLJ + fij * rxij;
                     fyijLJ = fyijLJ + fij * ryij;
                     fzijLJ = fzijLJ + fij * rzij;
             
                     if ( rij <= rCut )
                     {
                        fxijLJc = fxijLJc + fij * rxij;
                        fyijLJc = fyijLJc + fij * ryij;
                        fzijLJc = fzijLJc + fij * rzij;
                     }
             
                  }
/*           
 *                Coulomb forces   
 */          
                  if ( qi!=0.0 && qj!=0.0 )
                  {
                     fij  = fQQ * qi * qj / ( rij2 * rij );
                     fxijQQ = fxijQQ + fij * rxij;
                     fyijQQ = fyijQQ + fij * ryij;
                     fzijQQ = fzijQQ + fij * rzij;
                  }
/*           
 *             end of loops (i and j)   
 */          
               }
               }
/*           
 *             calculate the force along i->j direction   
 *             ffij > 0 :  repulsive force  
 *             ffij < 0 :  attractive force 
 */          
               ffijLJ = ( fxijLJ*cxij + fyijLJ*cyij + fzijLJ*czij ) / cij;
               ffijQQ = ( fxijQQ*cxij + fyijQQ*cyij + fzijQQ*czij ) / cij;
               ffijLJc = ( fxijLJc*cxij + fyijLJc*cyij + fzijLJc*czij ) / cij;
/*           
 *             determine rmin, rmax, bmin, bmax   
 */          
               cos_i = (rci2+cij2-rcj2)/(rci*cij*2.0);
               if ( fabs(cos_i)>1.0 ) 
               {
                  printf ("Error: cos_i<-1 or cos_i>1!\ncos_i = %f\n",cos_i);
                  printf ("rci = %f, cij = %f, rcj = %f\n",rci,cij,rcj);
                  exit(1);
               }
               h2 = rci2*(1-cos_i*cos_i);
               if ( fabs(rci2-rcj2) <= cij2 ) 
               {
                  rMin = sqrt(h2);
               } 
               else 
               {
                  rMin = rci<rcj?rci:rcj;
               }
               rMax = rci>rcj?rci:rcj;
               bMin = (int)( rMin / dR + 0.5 );
               bMax = (int)( rMax / dR + 0.5 );
               if ( bMin<0 ) 
               {
                  printf ("Error: bMin<0!\nrMin = %f, bMin = %d",rMin,bMin);
                  exit(1);
               }
               if ( bMax<0 ) 
               {
                  printf ("Error: bMax<0!\nrMax = %f, bMax = %d",rMax,bMax);
                  exit(1);
               }
               bMin = (bMin>=maxBin)?maxBin:bMin;
               bMax = (bMax>=maxBin)?maxBin:bMax;
/*           
 *             start of loop (bin)   
 */          
               for (bin=bMin; bin<bMax; bin++) 
               {
                  r  = ( (double)(bin) + 0.5 ) * dR;
                  r2 = r*r;
/*           
 *                determine number of intersections  
 */          
                  if ( (rci-r)*(rcj-r) < 0.0 ) 
                  {
                     inters = 1;
                  } 
                  else if ( rci>=r && rcj>=r ) 
                  {
                     if ( h2==r2 )
                     {
                        inters = 1;
                     }
                     else
                     {
                        inters = 2;
                     }
                  } 
                  else 
                  {
                     inters = 0;
                  }
                  cosine = sqrt(r2-h2)/r;
/*           
 *                add increment to pU   
 */          
                  if ( inters != 0 )
                  {
                     incrLJ = ffijLJ*cosine*inters/(4.0*pi*r2);
                     incrQQ = ffijQQ*cosine*inters/(4.0*pi*r2);
                     incrLJc = ffijLJc*cosine*inters/(4.0*pi*r2);
             
                     pULJ[bin] += incrLJ;
                     pUQQ[bin] += incrQQ;
                     pULJc[bin] += incrLJc;
             
                     averPULJ[bin] += incrLJ;
                     averPUQQ[bin] += incrQQ;
                     averPULJc[bin] += incrLJc;
                  }
/*           
 *             end of loop (bin)   
 */          
               }
/*
 *          end of loop (ires, jres)
 */
            }
            }
/*
 *       end of loop (pair)    
 */
         }
/*
 *       calculate work of formation  
 */
         workFull = 0.0;
         workCut = 0.0;
         for ( bin=0; bin<maxBin; bin++ ) 
         {
            r = ( (double)(bin) + 0.5 ) * dR;
            pK = kB*nA*temperature*averDens[bin];
            workFull += ( ( pK + pULJ[bin] + pUQQ[bin] ) * r*r ) * dR;
            workCut += ( ( pK + pULJc[bin] + pUQQ[bin] ) * r*r ) * dR;
         }
         workFull *= pi*2.0 ;
         workCut *= pi*2.0;
         printf ( "%8d%20.10f%20.10f\n", frame, workFull, workCut );
/*
 *    end of loop (frame)   
 */
      }
/*
 *    close trajectory file
 */
      xdrfile_close (xtc);

/*
 *    average arrays: averDens, averPU, molDens
 */
      for (bin=0; bin<maxBin; bin++) 
      {
         averDens[bin]  /= (nFrames-nStart);
         averPULJ[bin]  /= (nFrames-nStart);
         averPUQQ[bin]  /= (nFrames-nStart);
         averPULJc[bin] /= (nFrames-nStart);
         for ( mol=0; mol<molTypes; mol++ ) 
         {
            for ( res=0; res<resNr[mol]; res++ )
            {
               resDens[bin][mol][res] /= (nFrames-nStart);
            }
         }

      }

/*
 *    write results
 */
      printf ( "%s\n", "#Radius  P_N_full  P_N_cut  P_K  P_U_LJ  P_U_QQ  P_U_LJc" );
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         r = ( (double)(bin) + 0.5 ) * dR;
         pK = kB*nA*temperature*averDens[bin];
         printf ( "%8.3f%15.6f%15.6f%15.6f%15.6f%15.6f%15.6f\n", 
                  r, 
                  pK + averPULJ[bin] + averPUQQ[bin],
                  pK + averPULJc[bin] + averPUQQ[bin],
                  pK, 
                  averPULJ[bin],
                  averPUQQ[bin],
                  averPULJc[bin] );
      }
/*
 *    write molecular densities
 */
      printf ( "%s\n", "#Radius Total_Density Individual_Densities" );
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         r = ( (double)(bin) + 0.5 ) * dR;
         printf ( "%8.3f%20.10f", r, averDens[bin] );
         for ( mol=0; mol<molTypes; mol++ )
         {
            for ( res=0; res<resNr[mol]; res++ )
            {
               printf ( "%20.10f", resDens[bin][mol][res] );
            }
         }
         printf ( "\n" );
      }

/*
 *    free arrays
 */
      free(atom_in_mol);
      free(atom_in_res);
      free(molNr);

      free( charge  );
      free( mass    );
      free( sigma   );
      free( epsilon );

      free( iAtom );
      free( iMol  );

      free( molMass );

      free( mini );
      free( maxi );

      free( rx );
      free( ry );
      free( rz );

      free( cx );
      free( cy );
      free( cz );

      free( dens  );
      free( pULJ  );
      free( pUQQ  );
      free( pULJc );

      free( averDens  );
      free( averPULJ  );
      free( averPUQQ  );
      free( averPULJc );

      free( resDens );

/*****************************************************************************
 *    End of function main                                                   *
 *****************************************************************************/

   return 0;
   }

