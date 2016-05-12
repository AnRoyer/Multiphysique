# Solve a system using MUMPSseq in a sequential code

cmake_minimum_required(VERSION 2.7)

project(Multiphysique)

OPTION(MP_USE_MUMPS "use MUMPS solver instead of gmm::lu_solve"  OFF)

IF( (CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Intel") )
    ADD_DEFINITIONS(-std=gnu++0x)
ENDIF()

IF(MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_DEPRECATE)
ENDIF(MSVC)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")

include_directories(${PROJECT_SOURCE_DIR}/gmm)

IF(MP_USE_MUMPS)
    FIND_PACKAGE(MUMPS REQUIRED)
    include_directories(${MUMPS_INCLUDE_DIRS})

    IF(${BLA_VENDOR} MATCHES "OpenBlas")
    	find_library(LAPACK_LIBRARIES openblas)
    ELSE()
    	FIND_PACKAGE(LAPACK REQUIRED)   
    ENDIF()

    add_definitions(-DGMM_USES_LAPACK)  
    add_definitions(-DGMM_USES_MUMPS)
ENDIF()

# --------------------------------------------------------------
# Finds MPI (including MSMPI)
# --------------------------------------------------------------
if(WIN32)  
    FIND_PATH(MPI_INCLUDE_PATH NAMES  mpi.h HINTS "$ENV{MSMPI_INC}")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
        set(MS_MPI_ARCH_DIR x64)
        find_library(MPI_LIBRARIES msmpi PATHS "$ENV{MSMPI_LIB64}")
        message("64bits")
    else()
        set(MS_MPI_ARCH_DIR x86)
        find_library(MPI_LIBRARIES msmpi PATHS "$ENV{MSMPI_LIB32}")
        message("32bits")
    endif()
   IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        add_definitions(-DMSMPI_NO_SAL -D"__int64=long long") # mingw -sans visual
   ENDIF()
else()
    find_package(MPI REQUIRED)
endif()
# --------------------------------------------------------------
IF(WIN32)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_DEPRECATE)
ENDIF(WIN32)

include_directories(${MPI_INCLUDE_PATH})

FILE(GLOB SRCS *.h *.cpp)
add_executable(MP ${SRCS})

TARGET_LINK_LIBRARIES(MP ${MUMPS_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES})

