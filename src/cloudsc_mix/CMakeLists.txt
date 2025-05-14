# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# only double-precision currently supported
if( "${prec}" STREQUAL "dp" )

    # Define integer IDs corresponding to every language choice
    ## kernel language
    set(CUDA_LANG "1")
    set(HIP_LANG "2")
    set(SYCL_LANG "3")
    ## offload model
    set(ACC_OFFLOAD "1")
    set(OMP_OFFLOAD "2")

    # Select offload model
    if( NOT DEFINED CLOUDSC_GPU_OFFLOAD )
      set(CLOUDSC_GPU_OFFLOAD "ACC")
    endif()
    if( NOT ${CLOUDSC_GPU_OFFLOAD} MATCHES "ACC|OMP" )
      message(FATAL_ERROR "CLOUDSC_GPU_OFFLOAD: '${CLOUDSC_GPU_OFFLOAD}' is not a valid option! (Allowed: 'ACC' and 'OMP')")
    endif()

    set(GPU_OFFLOAD ${${CLOUDSC_GPU_OFFLOAD}_OFFLOAD})

    # Select kernel language
    if( NOT DEFINED CLOUDSC_GPU_LANG )
      set(CLOUDSC_GPU_LANG "CUDA")
    endif()
    if( NOT ${CLOUDSC_GPU_LANG} MATCHES "CUDA|HIP|SYCL" )
      message(FATAL_ERROR "CLOUDSC_GPU_LANG: '${CLOUDSC_GPU_LANG}' is not a valid option! (Allowed: 'CUDA', 'HIP' and 'SYCL')")
    endif()

    set(GPU_LANG ${${CLOUDSC_GPU_LANG}_LANG})

    if (CLOUDSC_GPU_LANG STREQUAL "CUDA")
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "-O3 -ffast-math")
    elseif (CLOUDSC_GPU_LANG STREQUAL "HIP")
	if(NOT DEFINED ROCM_PATH)
          if(DEFINED ENV{ROCM_PATH})
            set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCM has been installed")
          else()
            set(ROCM_PATH "/opt/rocm" CACHE PATH "Path to which ROCM has been installed")
          endif()
        endif()
	enable_language(HIP)
        find_package(hip REQUIRED)
	set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -O3 -ffast-math")
    elseif (CLOUDSC_GPU_LANG STREQUAL "SYCL")
	enable_language(CXX)
        set(CMAKE_CXX_STANDARD 17)
    endif()
    
    ecbuild_add_library(
        TARGET dwarf-cloudsc-gpu-lib-${prec}
        INSTALL_HEADERS LISTED
	TYPE STATIC
        SOURCES
    	    yoecldp_c.h
      	    cloudsc_c_k_caching.cpp
            cloudsc_c_k_caching_mod.F90
    	    PUBLIC_INCLUDES
    	    	$<INSTALL_INTERFACE:include>
            	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    	    PUBLIC_LIBS
                $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
		cloudsc-common-${prec}-lib
    )

    # Provide definitions to target
    foreach(_def CUDA_LANG HIP_LANG SYCL_LANG ACC_OFFLOAD OMP_OFFLOAD GPU_LANG GPU_OFFLOAD)
      target_compile_definitions(dwarf-cloudsc-gpu-lib-${prec} PUBLIC ${_def}=${${_def}})
    endforeach()

    target_include_directories(dwarf-cloudsc-gpu-lib-${prec} PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

    # CUDA specific
    if (CLOUDSC_GPU_LANG STREQUAL "CUDA")
      set_source_files_properties(cloudsc_c_k_caching.cpp PROPERTIES LANGUAGE CUDA)
      target_include_directories(
        dwarf-cloudsc-gpu-lib-${prec}
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
      )
      if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        target_compile_options(dwarf-cloudsc-gpu-lib-${prec} PRIVATE $<COMPILE_LANGUAGE:CUDA>)
      else()
        target_compile_options(dwarf-cloudsc-gpu-lib-${prec} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
	  -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
      endif()
      set_target_properties( dwarf-cloudsc-gpu-lib-${prec} PROPERTIES CUDA_SEPARABLE_COMPILATION ON) 
    # HIP specific
    elseif (CLOUDSC_GPU_LANG STREQUAL "HIP")
      set_source_files_properties(cloudsc_c_k_caching.cpp PROPERTIES LANGUAGE HIP)
      target_compile_options(dwarf-cloudsc-gpu-lib-${prec} PRIVATE --offload-arch=gfx90a) 
    # SYCL specific
    elseif (CLOUDSC_GPU_LANG STREQUAL "SYCL")
    
    else ()
    
    endif()
    
    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-mix-${prec}
	SOURCES dwarf_cloudsc_gpu.F90 cloudsc_driver_mod.F90
	LIBS dwarf-cloudsc-gpu-lib-${prec}
    )

    set_target_properties(dwarf-cloudsc-mix-${prec} PROPERTIES LINKER_LANGUAGE Fortran)

    # CUDA specific
    if (CLOUDSC_GPU_LANG STREQUAL "CUDA")
      set_target_properties( dwarf-cloudsc-mix-${prec} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    # HIP specific 
    elseif (CLOUDSC_GPU_LANG STREQUAL "HIP")
    # SYCL specific
    elseif (CLOUDSC_GPU_LANG STREQUAL "SYCL")
    else ()
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

endif()
