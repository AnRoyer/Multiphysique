# FindMUMPS.cmake - try to find MUMPS and associated libraries
#    (RoBo)
# ----------------------------------------------------------------------------   
# output:
#    MUMPS_FOUND             : TRUE/FALSE
#    MUMPS_LIBRARIES         : list of MUMPS libraries to link with
#    MUMPS_INCLUDE_DIRS      : where the *.h are                      [cached]
# ----------------------------------------------------------------------------


find_path(MUMPS_INCLUDE_DIRS NAMES "dmumps_c.h")

set(MUMPS_LIBRARIES "")
if(MSVC)
	SET(_MUMPS_LIBSET_SO dmumps_pord_seq_mkl)                        # mumps/seq (dynamic)
	SET(_MUMPS_LIBSET_A  dmumps_c dmumps_fortran 
                             mumps_common_c pord_c libseq_c 
                             libseq_fortran ifconsol libifcoremd 
                             libifportmd libmmd libirc svml_dispmd)  # mumps/seq (static)
else()
	SET(_MUMPS_LIBSET_SO dmumps_seq
                             dmumps_pord_seq_mkl 
                             dmumps_pord_seq_openblas)           # mumps/seq (dynamic)
	SET(_MUMPS_LIBSET_A  dmumps mumps_common pord mpiseq)    # mumps/seq (static)
endif()

# try the dyn libs first
find_library(FOUND_LIB_SO NAMES ${_MUMPS_LIBSET_SO})
#MESSAGE(STATUS "${_MUMPS_LIBSET_SO} = ${FOUND_LIB_SO}")
if(FOUND_LIB_SO)
	list(APPEND MUMPS_LIBRARIES ${FOUND_LIB_SO})
else()
	# ".so" not found => try ".a"
	foreach(LIB ${_MUMPS_LIBSET_A})
		find_library(FOUND_LIB_A ${LIB})
		MESSAGE(STATUS "${LIB}=${FOUND_LIB_A}")
		if(FOUND_LIB_A)
			list(APPEND MUMPS_LIBRARIES ${FOUND_LIB_A})
		endif()
		unset(FOUND_LIB_A CACHE)
	endforeach() 

	# add fortran libs
	IF(MUMPS_LIBRARIES)
            IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        	list(APPEND MUMPS_LIBRARIES -lgfortran)
	    ENDIF()
	    IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        	list(APPEND MUMPS_LIBRARIES -lifcoremt)
	    ENDIF()
        ENDIF()
endif() 
unset(FOUND_LIB_SO CACHE)


MESSAGE(STATUS "MUMPS_LIBRARIES = ${MUMPS_LIBRARIES}")

# ----------------------------------------------------------------------------   

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG 
				MUMPS_INCLUDE_DIRS
				MUMPS_LIBRARIES)

#MESSAGE(STATUS "MUMPS_FOUND = ${MUMPS_FOUND}")
# ----------------------------------------------------------------------------   
