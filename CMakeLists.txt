cmake_minimum_required (VERSION 3.15)

project (fit)
enable_language (CXX)
set (CMAKE_CXX_STANDARD 17)

set (WARNING_FLAGS "-Wall -Werror -Wextra -Wno-unused -Wno-sign-compare")
set (OPTIMIZATION_FLAGS "-O3")
set (DEBUG_FLAGS "-ggdb -Og")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_FLAGS} ${WARNING_FLAGS}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPTIMIZATION_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")

#set (MINUIT_STATIC_LIB "${CMAKE_SOURCE_DIR}/minuit/lib/libminuit.a")
add_subdirectory (minuit)
file (GLOB SRC ${CMAKE_SOURCE_DIR}/src/*.cpp)

include_directories (${CMAKE_SOURCE_DIR}/include)

add_library (obj_lib OBJECT ${SRC})

add_executable (fit ${CMAKE_SOURCE_DIR}/main.cpp $<TARGET_OBJECTS:obj_lib>)
target_link_libraries (fit -lgfortran minuit)

add_executable (checkmat ${CMAKE_SOURCE_DIR}/check-matrix.cpp $<TARGET_OBJECTS:obj_lib>)

install (TARGETS fit DESTINATION "${CMAKE_INSTALL_PREFIX}")
install (TARGETS checkmat DESTINATION "${CMAKE_INSTALL_PREFIX}")
