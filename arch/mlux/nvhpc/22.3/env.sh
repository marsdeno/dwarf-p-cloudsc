# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $1"
  module load $1
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

# Unload to be certain
module_unload Java
module_unload Python
module_unload HDF5
module_unload Boost
module_unload ParaStationMPI
module_unload NVHPC
module_unload CMake

# Load modules

module use /apps/USE/easybuild/release/2021.5/modules/all/

module purge
module load NVHPC/22.3

module_load psmpi/5.6.0-1-NVHPC-22.3
module_load CMake/3.20.4
module_load CUDA/11.3.1
module_load Boost/1.76.0-GCC-10.3.0

module_load Python/3.9.5-GCCcore-10.3.0-bare
module_load Java

export HDF5_ROOT=$HOME/dependencies/install

# Increase stack size to maximum
ulimit -S -s unlimited

# Restore tracing to stored setting
if [[ -n "$tracing_" ]]; then set -x; else set +x; fi

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
