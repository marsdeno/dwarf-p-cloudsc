# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS   "-fopenmp" CACHE STRING "" )
set( OpenMP_CXX_FLAGS   "-fopenmp" CACHE STRING "" )
set( OpenMP_Fortran_FLAGS   "-fopenmp" CACHE STRING "" )

####################################################################
# COMMON FLAGS
####################################################################

set(CMAKE_EXE_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")
