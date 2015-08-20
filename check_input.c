/*****************************************************************************
 This file is part of the SurfTensCalc program.                                      
 Copyright (C) 2012 Xin Li
                                                                           
 Filename:  check_input.c
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

/********************************************************************
 *    C program to check the input file                             *
 *    by Xin Li, TheoChem & Biology, KTH, Stockholm, Sweden         *
 *    Version 2012-03-16                                            *
 ********************************************************************/

   #include <math.h>
   #include <stdio.h>
   #include <stdlib.h>
   #include <xdrfile/xdrfile_trr.h>

/********************************************************************
 *    Main program                                                  *
 ********************************************************************/

   int main( int argc, char *argv[] )
   {

/********************************************************************
 *    Define variables                                              *
 ********************************************************************/

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
      double *averDens, *averPU, *dens, *pU ;
/*
 *    frame
 */
      int frame;
/*
 *    pair pointers for loop   
 */
      long int pair, totalPairs;
      int    *molI, *molJ;
/*
 *    center of mass   
 */
      double msum ;
/*
 *    variables for pair force calculation   
 */
      int    i, j, im, jm, iMin, iMax, jMin, jMax ;
      double mi, qi, sigi, epsi, qj, sigj, epsj;
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

      FILE    *file_par;
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
      char    line[256], tmp[256], filenameTRR[256];
      double  **mass, **charge, **sigma, **epsilon;
      int     k;
      int     nAtoms, nMols;
      double  *molMass;
/*
 *    open param.txt 
 */
      file_par = fopen("param.txt","r") ;
      if ( NULL==file_par ) {
         printf( "<> Error: Cannot open file param.txt !\n" ) ;
         exit(1);
      }
      mol = 0;
/*
 *    read comments (two lines)
 */
      fgets( line, sizeof( line ), file_par );
      fgets( line, sizeof( line ), file_par );
/*
 *    read name of trr file
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%s", tmp, filenameTRR );
/*
 *    read number of molecules
 */
      fgets( line, sizeof( line ), file_par );
      sscanf( line, "%s%d", tmp, &molTypes );
      atomNr = malloc(molTypes * sizeof(int));
      molNr  = malloc(molTypes * sizeof(int));

      printf ( "   There are %d types of molecules:\n", molTypes );
      printf ( "   %-16s", "MoleculeName" );

      for ( mol=0; mol<molTypes; mol++ ) 
      {
         fgets( line, sizeof( line ), file_par );
         sscanf( line, "%s%d", tmp, &molNr[mol] );

         printf ( "%8s", tmp );

      }

      printf ( "\n" );
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
/*
 *    print number of atoms in each type of molecule
 */
      printf ( "   %-16s", "ContainAtoms" );

      for ( mol=0; mol<molTypes; mol++ )
      {

         printf ( "%8d", atomNr[mol] );

      }

      printf ( "\n" );

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

      printf ( "   %-16s", "MolarMass" );

      for ( mol=0; mol<molTypes; mol++ ) 
      {
         molMass[mol] = 0.0;
         for ( atom=0; atom<atomNr[mol]; atom++ ) 
         {
            molMass[mol] += mass[mol][atom];
         }

         printf ( "%8.2f", molMass[mol] );

      }

      printf ( "\n" );

/*
 *    print number of atoms in each type of molecule
 */
      printf ( "   %-16s", "N_Molecules" );

      for ( mol=0; mol<molTypes; mol++ )
      {

         printf ( "%8d", molNr[mol] );

      }

      printf ( "\n" );

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

      x = calloc(nAtoms,sizeof(x[0]));
      v = calloc(nAtoms,sizeof(x[0]));
      f = calloc(nAtoms,sizeof(x[0]));

      totalPairs = nMols*(nMols-1)/2 ;
      molI = malloc( sizeof( int ) * totalPairs );
      molJ = malloc( sizeof( int ) * totalPairs );

/*
 *    read trr files for coordinates   
 */
      trr = xdrfile_open (filenameTRR,"r");
      read_trr_natoms (filenameTRR, &gmxNAtoms);
      /*printf("Total Number of Frames %10d\n", nFrames);*/
      /*printf("Starting from Frame    %10d\n", nStart);*/
/*
 *    start of loop (frame)   
 */
      for (frame=0; frame<2; frame++) 
      {
/*
 *       read trr file   
 */
         read_return = read_trr (trr,gmxNAtoms,&gmxStep,&gmxTime,&gmxLambda,gmxBox,x,v,f);
         if (read_return!=0) {
            printf ("Failed to read trr file!\n");
            exit(1);
         }
         if (gmxNAtoms!=nAtoms) {
            printf ("Incorrect number of Atoms!\n");
            exit(1);
         }
/*
 *       generate density profile   
 */
         for (im=0; im<nMols; im++) 
         {
            msum   = 0.0;
            iMin = mini[im];
            iMax = maxi[im];
            for (i=iMin; i<iMax; i++) 
            {
               mi = mass[iMol[i]][iAtom[i]];
               msum   += mi;
            }
            if ( fabs(msum-molMass[iMol[iMin]])>0.01 ) 
            {
               printf ("Incorrect molar weight of molecule %d!\n", im);
               exit(1);
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
 *    end of loop (frame)   
 */
      }

/*
 *    close trajectory file
 */
      xdrfile_close (trr);

      printf ( "   There are %d molecules in total.\n", nMols );
      printf ( "   There are %d atoms in total.\n", nAtoms );


/*****************************************************************************
 *    End of function main                                                   *
 *****************************************************************************/

   return 0;
   }

