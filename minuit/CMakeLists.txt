enable_language (Fortran)

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-second-underscore -malign-double -O3 -DCERNLIB_DOUBLE")
file (GLOB MINUIT_SRC ${CMAKE_SOURCE_DIR}/minuit/*.F)

add_library (minuit STATIC ${MINUIT_SRC})
target_include_directories (minuit PRIVATE ${CMAKE_SOURCE_DIR}/minuit)
