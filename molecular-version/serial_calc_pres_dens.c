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
      printf ( "%s\n", "#Frame Work_of_formation" );

/********************************************************************
 *    Define variables                                              *
 ********************************************************************/

/*
 *    number of frames in trajectory file 
 */
      int   frame, nFrames, nStart;
      nFrames = atoi(argv[1]);
      nStart = atoi(argv[2]);

      double temperature;
      temperature = atof(argv[3]); /* K */
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
 *    work of formation
 */
      double work;
/*
 *    coordinates of atoms and COMs, COM = center of mass 
 */
      double *rx, *ry, *rz;
      double *cx, *cy, *cz;
/*
 *    thickness of sperical layer, 0.3 Angstrom = 0.1 Sigma_Oxygen
 */  
      const double dR = 0.03 ;
      const int    maxBin = 1000 ;
      int    bin ;
      double dV ;
      double *averDens, *averPU, *dens, *pU ;
/*
 *    densities for different molecules 
 */
      double **molDens ;
/*
 *    increment, temporary variable 
 */
      double incr;
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
      double mi, qi, sigi, epsi, qj, sigj, epsj;
      double r, r2 ;
      double cxi, cyi, czi, rci2, rci ;
      double cxj, cyj, czj, rcj2, rcj ;
      double cxij, cyij, czij, cij2, cij ; 
      double rxij, ryij, rzij, rij2, rij ;
      double sr2, sr6, sr12, ffij, fxij, fyij, fzij, fij ;
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
/*
 *    atomNr = number of atoms in each type of molecule  
 *    molNr  = number of each type of molecule in the system 
 */
      int     *atomNr, *molNr, *iAtom, *iMol, *mini, *maxi;
      char    line[256], tmp[50], filenameXTC[50];
      double  **mass, **charge, **sigma, **epsilon;
      int     k;
      int     nAtoms, nMols;
      double  *molMass;
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
 *    read number of frames
 */
      fgets( line, sizeof( line ), file_par );
      /*sscanf( line, "%s%d", tmp, &nFrames );*/
      fgets( line, sizeof( line ), file_par );
      /*sscanf( line, "%s%d", tmp, &nStart );*/
      if ( nFrames<=nStart ) 
      {
         printf( "Error: nFrames is no larger than nStart!\n" ) ;
         exit(1);
      }
/*
 *    read name of xtc file
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%s", tmp, filenameXTC );
/*
 *    read number of molecules
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%d", tmp, &molTypes );
      atomNr = malloc(molTypes * sizeof(int));
      molNr  = malloc(molTypes * sizeof(int));
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &molNr[mol] );
      }
/*
 *    read atomic information
 */
      charge  = malloc(molTypes * sizeof(double *));
      mass    = malloc(molTypes * sizeof(double *));
      sigma   = malloc(molTypes * sizeof(double *));
      epsilon = malloc(molTypes * sizeof(double *));
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &atomNr[mol] );
         charge[mol]  = malloc(atomNr[mol] * sizeof(double));
         mass[mol]    = malloc(atomNr[mol] * sizeof(double));
         sigma[mol]   = malloc(atomNr[mol] * sizeof(double));
         epsilon[mol] = malloc(atomNr[mol] * sizeof(double));
         for ( atom=0; atom<atomNr[mol]; atom++ ) 
         {
            fgets( line, sizeof( line ), file_par );
            sscanf( line, "%lf%lf%lf%lf",
                    &charge[mol][atom], &mass[mol][atom],
                    &sigma[mol][atom],  &epsilon[mol][atom] );
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
         nAtoms += molNr[mol]*atomNr[mol];
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
         for ( atom=0; atom<atomNr[mol]; atom++ ) 
         {
            molMass[mol] += mass[mol][atom];
         }
      }
/*
 *    assign atom index and mol index for each atom in the system 
 *    i = atom index, j = mol index
 */
      i = 0;
      j = 0;
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         for ( k=0; k<molNr[mol]; k++ ) 
         {
            for ( atom=0; atom<atomNr[mol]; atom++ ) 
            {
               iMol[i] = mol;
               iAtom[i] = atom;
               i++;
            }
            j++;
         }
      }

      if ( nAtoms != i ) 
      {
        printf("Incorrect number of atoms in the system!\n");
        printf("Program stops at assigning atom index and mol index!\n");
        exit(1);
      }

      if ( nMols != j ) 
      {
        printf("Incorrect number of molecules in the system!\n");
        printf("Program stops at assigning atom index and mol index!\n");
        exit(1);
      }
/*
 *    assign starting index and ending index for each molecule in the system 
 *    i = atom index, j = mol index
 */
      mini = malloc(nMols * sizeof(int));
      maxi = malloc(nMols * sizeof(int));
      i = 0;
      j = 0;
      for ( mol=0; mol<molTypes; mol++ ) 
      {
         for ( k=0; k<molNr[mol]; k++ ) 
         {
            mini[j] = i;
            i += atomNr[mol];
            maxi[j] = i;
            j++;
         }
      }

      if ( nAtoms != i ) 
      {
        printf("Incorrect number of atoms in the system!\n");
        printf("Program stops at assigning atom index and mol index!\n");
        exit(1);
      }

      if ( nMols != j ) 
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

      cx = malloc( sizeof(double) * nMols );
      cy = malloc( sizeof(double) * nMols );
      cz = malloc( sizeof(double) * nMols );
      if ( NULL==cx || NULL==cy || NULL==cz ) 
      {
         printf ("Unable to allocate arrays: cx, cy, cz!\n");
         exit(1);
      }

      averDens  = malloc( sizeof( double ) * maxBin );
      averPU    = malloc( sizeof( double ) * maxBin );
      if ( NULL==averDens || NULL==averPU ) 
      {
         printf ("Unable to allocate arrays: averDens, averPU!\n");
         exit(1);
      }

      dens  = malloc( sizeof( double ) * maxBin );
      pU    = malloc( sizeof( double ) * maxBin );
      if ( NULL==dens || NULL==pU ) 
      {
         printf ("Unable to allocate arrays: dens, pU!\n");
         exit(1);
      }

      molDens = malloc(maxBin * sizeof(double *));
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         molDens[bin] = malloc(molTypes * sizeof(double));
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
         averPU[bin]  = 0.0;
         for ( mol=0; mol<molTypes; mol++ ) 
         {
            molDens[bin][mol] = 0.0;
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
            pU[bin] = 0.0;
         }
/*
 *       generate density profile   
 */
         for (im=0; im<nMols; im++) 
         {
            cx[im] = 0.0;
            cy[im] = 0.0;
            cz[im] = 0.0;
            msum   = 0.0;
            iMin = mini[im];
            iMax = maxi[im];
            for (i=iMin; i<iMax; i++) 
            {
               mi = mass[iMol[i]][iAtom[i]];
               cx[im] += mi*rx[i];
               cy[im] += mi*ry[i];
               cz[im] += mi*rz[i];
               msum   += mi;
            }
            if ( fabs(msum-molMass[iMol[iMin]])>0.01 ) 
            {
               printf ("Incorrect molar weight of molecule %d!\n", im);
               exit(1);
            }
            cx[im] /= msum;
            cy[im] /= msum;
            cz[im] /= msum;
            rci2 = cx[im]*cx[im] + cy[im]*cy[im] + cz[im]*cz[im];
            rci = sqrt( rci2 );
            bin = (int)( rci/dR );
            if ( bin < maxBin ) 
            {
/*
 *             volume of the spherical layer 
 */
               dV = pi*4.0/3.0*( pow((double)(bin+1)*dR,3) - pow((double)bin*dR,3) );
/*
 *             add increment to dens[bin] and averDens[bin]
 */
               incr = 1.0/dV;
               dens[bin] += incr;
               averDens[bin] += incr;
               molDens[bin][iMol[iMin]] += incr;
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
/*
 *          center of mass distances   
 */
            rci2 = cx[im]*cx[im] + cy[im]*cy[im] + cz[im]*cz[im];
            rci = sqrt( rci2 );

            rcj2 = cx[jm]*cx[jm] + cy[jm]*cy[jm] + cz[jm]*cz[jm];
            rcj = sqrt( rcj2 );

            cxij = cx[jm] - cx[im];
            cyij = cy[jm] - cy[im];
            czij = cz[jm] - cz[im];
            cij2 = cxij*cxij + cyij*cyij + czij*czij;
            cij  = sqrt( cij2 );
/*
 *          calculate forces between iw and jw  
 */
            fxij = 0.0;
            fyij = 0.0;
            fzij = 0.0;
/*
 *          determine the upper and lower boundary of i and j   
 */
            iMin = mini[im];
            iMax = maxi[im];
            jMin = mini[jm];
            jMax = maxi[jm];
/*
 *          start of loops (i and j)   
 */
            for (i=iMin; i<iMax; i++) 
            {
            for (j=jMin; j<jMax; j++) 
            {
/*
 *             interatomic distances   
 */
               rxij = rx[j] - rx[i];
               ryij = ry[j] - ry[i];
               rzij = rz[j] - rz[i];
               rij2 = rxij*rxij + ryij*ryij + rzij*rzij;
               rij  = sqrt( rij2 );
/*
 *             parameters for atom i and j
 */
               mol  = iMol[i];
               atom = iAtom[i];
               qi   =  charge[mol][atom];
               sigi =   sigma[mol][atom];
               epsi = epsilon[mol][atom];
               mol  = iMol[j];
               atom = iAtom[j];
               qj   =  charge[mol][atom];
               sigj =   sigma[mol][atom];
               epsj = epsilon[mol][atom];
/*
 *             vdw forces, LJ potential ( sigma, epsilon )
 *             U = 4 * eps * ( sig^12/R^12 - sig^6/R^6 )
 *             F = 4 * eps * ( 12*sig^12/R^13 - 6*sig^6/R^7 )
 *             vecF = 4 * eps * ( 12*sig^12/R^14 - 6*sig^6/R^8 ) * vecR
 *                  = 24 * eps * ( 2*(sig/R)^12 - (sig/R)^6 ) / R^2 * vecR
 *
 *             OPLS combination rule:
 *             sig(ij) = sqrt( sig(i) * sig(j) )
 *             eps(ij) = sqrt( eps(i) * eps(j) )
 */
               if ( epsi>0.0 && epsj>0.0 ) 
               {
                  sr2  = sigi*sigj/rij2;
                  sr6  = sr2 * sr2 * sr2;
                  sr12 = sr6 * sr6;
                  fij  = (sr12*2.0-sr6)/rij2;
                  fij  = fij * sqrt(epsi*epsj) * 24.0;
                  fxij = fxij + fij * rxij;
                  fyij = fyij + fij * ryij;
                  fzij = fzij + fij * rzij;
               }
/*
 *             Coulomb forces   
 */
               fij  = fQQ * qi * qj / ( rij2 * rij );
               fxij = fxij + fij * rxij;
               fyij = fyij + fij * ryij;
               fzij = fzij + fij * rzij;
/*
 *          end of loops (i and j)   
 */
            }
            }
/*
 *          calculate the force along i->j direction   
 *          ffij > 0 :  repulsive force  
 *          ffij < 0 :  attractive force 
 */
            ffij = ( fxij*cxij + fyij*cyij + fzij*czij ) / cij;

/*
 *          determine rmin, rmax, bmin, bmax   
 */
            cos_i = (rci2+cij2-rcj2)/(rci*cij*2.0);
            if ( cos_i<-1 || cos_i>1 ) 
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
 *          start of loop (bin)   
 */
            for (bin=bMin; bin<bMax; bin++) 
            {
               r  = ( (double)(bin) + 0.5 ) * dR;
               r2 = r*r;
/*
 *             determine number of intersections  
 */
               if ( (rci-r)*(rcj-r) < 0.0 ) 
               {
                  inters = 1;
               } 
               else if ( rci>=r && rcj>=r ) 
               {
                  inters = 2;
               } 
               else 
               {
                  inters = 0;
               }
               cosine = sqrt(r2-h2)/r;
/*
 *             add increment to pU   
 */
               if ( inters != 0 )
               {
                  incr = ffij*cosine*inters/(4.0*pi*r2);
                  pU[bin] += incr;
                  averPU[bin] += incr;
               }
/*
 *          end of loop (bin)   
 */
            }
/*
 *       end of loop (pair)    
 */
         }
/*
 *       calculate work of formation  
 */
         work = 0.0;
         for ( bin=0; bin<maxBin; bin++ ) 
         {
            r = ( (double)(bin) + 0.5 ) * dR;
            work += ( ( kB*nA*temperature*dens[bin] + pU[bin] ) * r*r ) * dR;
         }
         work *= pi*2.0 ;
         printf ( "%8d%20.10f\n", frame, work );
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
         averPU[bin]    /= (nFrames-nStart);
         for ( mol=0; mol<molTypes; mol++ ) 
         {
            molDens[bin][mol] /= (nFrames-nStart);
         }

      }

/*
 *    write results
 */
      printf ( "%s\n", "#Radius Pressure_K Pressure_U Pressure_N" );
      for ( bin=0; bin<maxBin; bin++ ) 
      {
         r = ( (double)(bin) + 0.5 ) * dR;
         printf ( "%8.3f%20.10f%20.10f%20.10f\n", r, 
                  kB*nA*temperature*averDens[bin], averPU[bin],
                  kB*nA*temperature*averDens[bin] + averPU[bin] );
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
            printf ( "%20.10f", molDens[bin][mol] );
         }
         printf ( "\n" );
      }

/*
 *    free arrays
 */
      free(atomNr);
      free(molNr );

      for ( mol=0; mol<molTypes; mol++ ) 
      {
        free( charge[mol]  );
        free( mass[mol]    );
        free( sigma[mol]   );
        free( epsilon[mol] );
      }
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

      free( averDens  );
      free( averPU  );

      for ( bin=0; bin<maxBin; bin++ ) 
      {
        free( molDens[bin] );
      }
      free( molDens );

/*****************************************************************************
 *    End of function main                                                   *
 *****************************************************************************/

   return 0;
   }

