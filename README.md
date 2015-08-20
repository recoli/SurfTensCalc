# SurfTensCalc
An old code that calculates surface tension of a spherical water droplet.

Requires the XTC Library
http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library

Usage:

mpirun -np 40 ./mpi-main.x 100001 50001 298.0 2.5 10.88960
./sphere_process_data.x 40 10.88960 > surftens.dat

Note:
a) 40 is number of CPUs
b) 100001 50001 means the calculation starts from frame 50001 and ends at frame 100000
c) 298.0 is temperature (K)
d) 2.5 is cutoff radius (nm)
e) 10.88960 is length of the simulation box (nm)
