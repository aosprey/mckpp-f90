# MC-KPP
 
Multi-Column configuration of the K Profile Parameterisation mixed-layer ocean model. 

This is a development version. We are working towards a modern f90 MPI parallel version of the original: http://puma.nerc.ac.uk/trac/KPP_ocean

See the wiki for more information on the project: https://github.com/aosprey/mckpp-f90/wiki

## Installation ## 

Pre-requisites: 
* Fortran 90 compiler
* Netcdf 
* XIOS parallel I/O library: https://forge.ipsl.jussieu.fr/ioserver
* FCM build system: http://metomi.github.io/fcm/doc/

Compiling: 
* Create an FCM configuration file with build flags for your platform. Use one of the existing ones as an example. 
* Build with fcm make -f <new-config-file> 
  
Running: 
* 
