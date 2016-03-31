# Solve a system using MUMPSseq in a sequential code

cmake_minimum_required(VERSION 2.7)

project(Multiphysique)

OPTION(MP_USE_MUMPS "use MUMPS solver instead of gmm::lu_solve"  OFF)

ADD_DEFINITIONS(-std=gnu++0x)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")

include_directories(${PROJECT_SOURCE_DIR}/gmm)

IF(MP_USE_MUMPS)
    FIND_PACKAGE(MUMPS REQUIRED)
    include_directories(${MUMPS_INCLUDE_DIRS})

    FIND_PACKAGE(LAPACK REQUIRED)

    add_definitions(-DGMM_USES_LAPACK)  
    add_definitions(-DGMM_USES_MUMPS)
ENDIF()

FILE(GLOB SRCS *.h *.cpp)
add_executable(MP ${SRCS})

TARGET_LINK_LIBRARIES(MP ${MUMPS_LIBRARIES} ${LAPACK_LIBRARIES})

