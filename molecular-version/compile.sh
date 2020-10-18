#!/bin/bash

# compile programs
  mpicc -O2 -o mpi-main.x mpi-main.c
  gcc -O2 -o check_input.x check_input.c -I./xdrfile/include -L./xdrfile/lib -lxdrfile
  gcc -O2 -o serial_calc_pres_dens.x serial_calc_pres_dens.c -I./xdrfile/include -L./xdrfile/lib -lxdrfile -lm
  gcc -O2 -o serial_process_data.x serial_process_data.c -lm
  gcc -O2 -o fit_dens.x fit_dens.c -lm
