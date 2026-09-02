cmake_minimum_required(VERSION 3.22)
project(jlinalg_cholmod LANGUAGES C)

find_package(JNI REQUIRED)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SUITESPARSE_ENABLE_PROJECTS cholmod CACHE STRING "" FORCE)
set(SUITESPARSE_USE_CUDA OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_USE_OPENMP OFF CACHE BOOL "" FORCE)
set(CHOLMOD_USE_CUDA OFF CACHE BOOL "" FORCE)
set(CHOLMOD_USE_OPENMP OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_DEMOS OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

if(DEFINED JLINALG_MKL_ROOT)
    file(TO_CMAKE_PATH "${JLINALG_MKL_ROOT}" JLINALG_MKL_ROOT_CMAKE)
    set(BLA_VENDOR Intel10_64lp CACHE STRING "" FORCE)
    set(BLAS_LIBRARIES
        "${JLINALG_MKL_ROOT_CMAKE}/lib/intel64/mkl_rt.lib"
        CACHE FILEPATH "" FORCE)
    set(LAPACK_LIBRARIES "${BLAS_LIBRARIES}" CACHE FILEPATH "" FORCE)
    set(BLAS_INCLUDE_DIRS "${JLINALG_MKL_ROOT_CMAKE}/include"
        CACHE PATH "" FORCE)
    if(MSVC)
        string(APPEND CMAKE_C_FLAGS
            " /I\"${JLINALG_MKL_ROOT_CMAKE}/include\"")
    endif()
endif()

include(FetchContent)
FetchContent_Declare(SuiteSparse
    URL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.12.2.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    URL_HASH SHA256=679412daa5f69af96d6976595c1ac64f252287a56e98cc4a8155d09cc7fd69e8)
FetchContent_MakeAvailable(SuiteSparse)

add_library(jlinalg_cholmod SHARED cholmod_jni.c)
target_include_directories(jlinalg_cholmod PRIVATE ${JNI_INCLUDE_DIRS})
target_link_libraries(jlinalg_cholmod PRIVATE CHOLMOD_static ${JNI_LIBRARIES})
target_compile_definitions(jlinalg_cholmod PRIVATE
    JLINALG_CHOLMOD_VERSION="SuiteSparse-7.12.2/CHOLMOD-5.3.4")
set_target_properties(jlinalg_cholmod PROPERTIES
    C_STANDARD 11
    C_STANDARD_REQUIRED ON
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
