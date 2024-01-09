# MC-KPP
 
Multi-Column configuration of the K Profile Parameterisation mixed-layer ocean model. 

This is a development version. We are working towards a modern Fortran, MPI parallel version of the original: https://github.com/NCAS-CMS/KPP_ocean

See the wiki for more information on the project: https://github.com/aosprey/mckpp-f90/wiki

## Installation ## 

Pre-requisites: 
* Modern Fortran compiler (tested gfortran v11.4.0 and Cray fortran v15.0.0) 
* XIOS parallel I/O library: https://forge.ipsl.jussieu.fr/ioserver (tested r2245)
  * This needs netcdf4 and hdf5
* FCM build system: http://metomi.github.io/fcm/doc/

Compiling: 
* Create an FCM build config file for your compiler and platform, using one of the existing ones as template. 
* For ARCHER2 use `fcm-make-crayftn-archer2.cfg`, and load the following modules: 
  ```
  module load cray-hdf5-parallel
  module load cray-netcdf-hdf5parallel
  ```
* Then build with `fcm make -f <config-file>`
  
Instructions for running a forced test on ARCHER2: 
* Create a run directory. 
* Copy the control scripts and the executable into the run directory: 
  ```
  cp ../mckpp-f90/run/* . 
  cp ../mckpp-f90/build/bin/KPP_ocean . 
  ```
* Sym-link the input files to the run directory
  ```
  ln -s $UMDIR/kpp/terramaris_forced/*.nc . 
  ```
* Check the Slurm batch settings and submit: 
  ```
  sbatch KPPocean.slurm
  ```
