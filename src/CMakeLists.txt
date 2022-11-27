cmake_minimum_required(VERSION 3.21)

project(hewes LANGUAGES C CXX)
set(PROJECT_VERSION "0.2.2")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: \
Debug Release RelWithDebInfo MinSizeRel." FORCE)
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Find MPI and optionally OpenMP
find_package(MPI 3.1 REQUIRED COMPONENTS CXX)
find_package(OpenMP 4.5 COMPONENTS C CXX)

# Find SUNDIALS
# Set path via `cmake -D SUNDIALS_DIR=/path/to/sundials/installation` or here
if(NOT DEFINED SUNDIALS_DIR)
    set(SUNDIALS_DIR "/Users/andi/Documents/sundials-6.3.0/INSTDIR_Rel")
    #set(SUNDIALS_DIR "/dss/dsshome1/lxc0A/ru68dab/Downloads/sundials-6.2.0/INSTDIR_Rel")
endif()
# Search only in given path to avoid unintended package manager installation
find_package(SUNDIALS 6.0 NO_DEFAULT_PATH PATHS ${SUNDIALS_DIR})
if(${SUNDIALS_FOUND})
    message(STATUS "Found SUNDIALS in ${SUNDIALS_DIR}.")
else()
    message(WARNING "SUNDIALS not found. Fetching CVODE ...")
    include(FetchContent)
    set(BUILD_STATIC_LIBS OFF CACHE INTERNAL "")
    set(ENABLE_MPI ON CACHE INTERNAL "")
    set(ENABLE_OPENMP ON CACHE INTERNAL "")
    set(EXAMPLES_ENABLE_C OFF CACHE INTERNAL "")
    set(EXAMPLES_INSTALL OFF CACHE INTERNAL "")
    set(BUILD_ARKODE OFF CACHE INTERNAL "")
    set(BUILD_CVODES OFF CACHE INTERNAL "")
    set(BUILD_IDA OFF CACHE INTERNAL "")
    set(BUILD_IDAS OFF CACHE INTERNAL "")
    set(BUILD_KINSOL OFF CACHE INTERNAL "")
    FetchContent_Declare(SUNDIALS
            GIT_REPOSITORY      https://github.com/LLNL/sundials
            GIT_TAG             v6.3.0)
    FetchContent_MakeAvailable(SUNDIALS)
endif()

# Project sources
list(APPEND sources
    main.cpp
    DerivationStencils.cpp
    ICSetters.cpp
    LatticePatch.cpp
    Outputters.cpp
    SimulationClass.cpp
    SimulationFunctions.cpp
    TimeEvolutionFunctions.cpp)

# Simulation executable
add_executable(hewes ${sources})

target_compile_options(hewes PUBLIC -O2)

# Link Sundials, MPI, and OpenMP publicly
target_link_libraries(hewes PRIVATE SUNDIALS::cvode
    SUNDIALS::nvecopenmp SUNDIALS::nvecmpiplusx
    SUNDIALS::nvecparallel MPI::MPI_CXX)

if(${OpenMP_FOUND})
    target_link_libraries(hewes PRIVATE OpenMP::OpenMP_CXX)
endif()

