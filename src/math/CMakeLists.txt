set( FLAGS_DEBUG -g -fopenmp -Wall -Wextra -Wno-unused-parameter )
set( FLAGS_RELEASE -O3 -fopenmp -Wall -Wextra -Wno-unused-parameter )

set( LIB_SOURCES Mat.cpp Sparse.cpp Vector.cpp )
set( LIB_HEADERS Mat.h Sparse.h Vector.h )

add_library(
    smath
    ${LIB_SOURCES}
    ${LIB_HEADERS}
)

target_compile_options(
    smath
    PUBLIC
        ${FLAGS_DEBUG}
)

set_target_properties(
        smath
        PROPERTIES CXX_STANDARD 11
)
