cmake_minimum_required(VERSION 3.0)

# set the project name
project(MO815_tasks VERSION 0.1)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_VERBOSE_MAKEFILE "ON")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=gnu11 -Wall -Wno-unused-result -pedantic")

file(GLOB EXTERNAL
        ${CMAKE_CURRENT_LIST_DIR}/externals/*/)

foreach(lib ${EXTERNAL})
    if(EXISTS ${lib}/src)
        message(STATUS ${lib};)
        list(APPEND EXTERNAL_INCLUDE ${lib}/include)
        list(APPEND EXTERNAL_SRC ${lib}/src)
        list(APPEND EXTERNAL_LIB ${lib}/lib)
        message(STATUS ${EXTERNAL_SRC};)
    endif()
endforeach()

foreach(ext_src ${EXTERNAL_SRC})
    aux_source_directory(${ext_src} IFT_SRC)
endforeach()

foreach(ext_inc ${EXTERNAL_INCLUDE})
    include_directories(${ext_inc})
endforeach()

#add libraries

add_library(ift STATIC IMPORTED)
set_target_properties(ift PROPERTIES
    IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/libift.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include"
)

include_directories($ENV{NEWIFT_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS} -fsanitize=address")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package( BLAS REQUIRED )
include_directories(${BLAS_INCLUDE_DIR})
find_package(LAPACK)

# add the executable
add_executable(getslice getslice.c)
target_link_libraries(getslice ift)
target_link_libraries(getslice m)
target_link_libraries(getslice blas ${BLAS_LIBRARIES})
target_link_libraries(getslice -lstdc++)
target_link_libraries(getslice ${LAPACK_LIBRARIES})

# add the executable
add_executable(MIP MIP.c)
target_link_libraries(MIP ift)
target_link_libraries(MIP m)
target_link_libraries(MIP blas ${BLAS_LIBRARIES})
target_link_libraries(MIP -lstdc++)
target_link_libraries(MIP ${LAPACK_LIBRARIES})

# add the executable
add_executable(reslicing reslicing.c)
target_link_libraries(reslicing ift)
target_link_libraries(reslicing m)
target_link_libraries(reslicing blas ${BLAS_LIBRARIES})
target_link_libraries(reslicing -lstdc++)
target_link_libraries(reslicing ${LAPACK_LIBRARIES})

# add the executable
add_executable(watershed watershed.c)
target_link_libraries(watershed ift)
target_link_libraries(watershed m)
target_link_libraries(watershed blas ${BLAS_LIBRARIES})
target_link_libraries(watershed -lstdc++)
target_link_libraries(watershed ${LAPACK_LIBRARIES})

# add the executable
add_executable(o_watershed orientedwatershed.c)
target_link_libraries(o_watershed ift)
target_link_libraries(o_watershed m)
target_link_libraries(o_watershed blas ${BLAS_LIBRARIES})
target_link_libraries(o_watershed -lstdc++)
target_link_libraries(o_watershed ${LAPACK_LIBRARIES})