# Sequential Gromacs
#-------------------
# Modules:
# cmake

cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/home/ovillegas/.applications/gromacs

# Gromacs in Parallel
#--------------------
# Modules:
# MPI 2.0 standard
# cmake

cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/home/ovillegas/.applications/gromacs_mpi -DGMX_MPI=on
