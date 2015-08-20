#!/bin/bash -x

# compile programs
  mpicc -O2 -o mpi-main.x mpi-main.c 
  gcc -O2 -o check_input.x check_input.c -lxdrfile
  gcc -O2 -o sphere_calc_pres_dens.x sphere_calc_pres_dens.c -lxdrfile
  gcc -O2 -o sphere_process_data.x sphere_process_data.c

