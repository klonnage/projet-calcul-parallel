add_subdirectory(math)

set(PROJ_SOURCES functions.cpp GradConj.cpp Input.cpp main.cpp)
set(PROJ_HEADERS functions.h Input.h GradConj.h )

add_executable(
    exe
    ${PROJ_SOURCES}
    ${PROJ_HEADERS}
)

SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -L$ENV{MKLROOT}/lib/intel64 -Wl --no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl" )

target_include_directories(
    exe
    PRIVATE
        "${PROJECT_SOURCE_DIR}/src/math"
)

target_link_libraries(
    exe
    PUBLIC 
        smath
)

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