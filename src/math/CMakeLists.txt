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
         -DMKL_ILP64 -m64 -I${MKLROOT}/include
)

add_link_options( -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl)

set_target_properties(
        smath
        PROPERTIES CXX_STANDARD 11
)
