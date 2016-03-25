# Solve a system using MUMPSseq in a sequential code

cmake_minimum_required(VERSION 2.7)

project(Multiphysique)

ADD_DEFINITIONS(-std=gnu++0x)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
FIND_PACKAGE(MUMPS REQUIRED)
include_directories(${MUMPS_INCLUDE_DIRS})

FIND_PACKAGE(LAPACK REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/gmm)
add_definitions(-DGMM_USES_LAPACK)  
add_definitions(-DGMM_USES_MUMPS)

FILE(GLOB SRCS *.h *.cpp)
add_executable(MP ${SRCS})

TARGET_LINK_LIBRARIES(MP ${MUMPS_LIBRARIES} ${LAPACK_LIBRARIES})

