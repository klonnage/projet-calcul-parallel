#add_subdirectory(math)

set(PROJ_SOURCES functions.cpp GradConj.cpp Input.cpp main.cpp math/Sparse.cpp math/Vector.cpp math/Mat.cpp)
set(PROJ_HEADERS functions.h Input.h GradConj.h math/Sparse.h math/Vector.h math/Mat.h)

set( FLAGS_RELEASE -O3 -fopenmp -Wall -Wextra -Wno-unused-parameter )

add_executable(
    exe
    ${PROJ_SOURCES}
    ${PROJ_HEADERS}
)


target_include_directories(
    exe
    PRIVATE
        "${PROJECT_SOURCE_DIR}/src/math"
)

#target_link_libraries(
#    exe
#    PUBLIC 
#        smath
#)

target_compile_options(
    exe
    PUBLIC
        ${FLAGS_DEBUG}
)

set_target_properties(
        exe
        PROPERTIES CXX_STANDARD 11
)

find_package(MPI)
    if (MPI_CXX_FOUND)
        target_link_libraries(exe PUBLIC MPI::MPI_CXX ${MPI_CXX_LINK_FLAGS})
    else()
        message(FATAL_ERROR "MPI not found !")
endif()

target_link_options(exe PUBLIC -L$ENV{MKLROOT}/lib/intel64 -Wl --no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl )
