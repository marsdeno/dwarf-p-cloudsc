# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

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

module_load llvm/EPI-development
#module_load openmpi/riscv/4.1.6_llvm1.0
module_load openmpi/ubuntu/4.1.6_llvm1.0
module_load hdf5/ubuntu/1.14.6_llvmEPIdev
module_load cmake/3.28.1

export FC=flang
export CC=clang
export CXX=clang++
set -x

export HDF5_PATH=$HDF5_INSTALL_DIR
export MPI_HOME=/apps/riscv/ubuntu/openmpi/4.1.6_llvm1.0/
export CPATH=/apps/riscv/ubuntu/openmpi/4.1.6_llvm1.0/include/:$CPATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
