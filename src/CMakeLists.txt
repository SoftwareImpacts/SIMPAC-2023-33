cmake_minimum_required(VERSION 3.21)

project(HEWES LANGUAGES C CXX VERSION 0.2.3)

# Print infos
include(CMakePrintHelpers)

cmake_print_variables(CMAKE_PROJECT_NAME)
cmake_print_variables(CMAKE_PROJECT_VERSION)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_print_variables(CMAKE_CXX_COMPILER_ID)
cmake_print_variables(CMAKE_CXX_COMPILER)
cmake_print_variables(CMAKE_CXX_STANDARD)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: \
Debug Release RelWithDebInfo MinSizeRel." FORCE)
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

option(BUILD_SHARED_LIBS "Build shared libraries." ON)

# Use of MPI
find_package(MPI 3.1 COMPONENTS CXX)
if(${MPI_FOUND} AND NOT DEFINED USE_MPI)
    set(USE_MPI ON)
    message(STATUS "Using MPI.")
    cmake_print_variables(MPI_INCLUDE_PATH)
    cmake_print_variables(MPI_CXX_LIBRARIES)
endif()

# Use of OpenMP
find_package(OpenMP 4.5 COMPONENTS C CXX)  # C required for Sundials
if(${OpenMP_FOUND} AND NOT DEFINED USE_OPENMP)
    set(USE_OPENMP ON)
    message(STATUS "Using OpenMP.")
    cmake_print_variables(OpenMP_CXX_INCLUDE_DIRS)
    cmake_print_variables(OpenMP_CXX_LIBRARIES)
endif()

# Find SUNDIALS
# Set path via `cmake -D SUNDIALS_DIR=/path/to/sundials/installation` or here
if(NOT DEFINED SUNDIALS_DIR)
    set(SUNDIALS_DIR "~/Documents/sundials-6.3.0/INSTDIR_Rel")
    #set(SUNDIALS_DIR "~/Downloads/sundials-6.2.0/INSTDIR_Rel")
endif()
# Search only in given path to avoid unintended package manager installation
find_package(SUNDIALS 6.0 NO_DEFAULT_PATH PATHS ${SUNDIALS_DIR})
if(${SUNDIALS_FOUND})
    message(STATUS "Found SUNDIALS in ${SUNDIALS_DIR}.")
else()
    message(WARNING "SUNDIALS not found. Fetching CVODE ...")
    include(FetchContent)
    set(BUILD_STATIC_LIBS OFF CACHE INTERNAL "")
    if(${MPI_FOUND})
        set(ENABLE_MPI ON CACHE INTERNAL "")
    endif()
    if(${OpenMP_FOUND})
        set(ENABLE_OPENMP ON CACHE INTERNAL "")
    endif()
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

# Exclude MPI stuff if not used
if(${MPI_FOUND} AND ${USE_MPI})
    target_compile_definitions(hewes PUBLIC "_MPI")
endif()

target_compile_options(hewes PUBLIC -O2)

# Link Sundials, MPI, and OpenMP
target_link_libraries(hewes PRIVATE
    SUNDIALS::cvode)

if(${MPI_FOUND} AND ${USE_MPI})
    target_link_libraries(hewes PRIVATE
        MPI::MPI_CXX
        SUNDIALS::nvecparallel)
endif()

if(${OpenMP_FOUND} AND ${USE_OPENMP})
    target_link_libraries(hewes PRIVATE
        OpenMP::OpenMP_CXX
        SUNDIALS::nvecopenmp)
endif()

if(${MPI_FOUND} AND ${OpenMP_FOUND} AND ${USE_MPI} AND ${USE_OPENMP})
    target_link_libraries(hewes PRIVATE
        SUNDIALS::nvecmpiplusx)
endif()

# Print system info
message(STATUS "\nSystem info:")
foreach(key
  IN ITEMS
    OS_NAME
    OS_RELEASE
    OS_VERSION
    OS_PLATFORM
    NUMBER_OF_LOGICAL_CORES
    NUMBER_OF_PHYSICAL_CORES
    TOTAL_VIRTUAL_MEMORY
    AVAILABLE_VIRTUAL_MEMORY
    TOTAL_PHYSICAL_MEMORY
    AVAILABLE_PHYSICAL_MEMORY
    IS_64BIT
    HAS_FPU
    HAS_MMX
    HAS_MMX_PLUS
    HAS_SSE
    HAS_SSE2
    HAS_SSE_FP
    HAS_SSE_MMX
    HAS_AMD_3DNOW
    HAS_AMD_3DNOW_PLUS
    HAS_IA64
  )
  cmake_host_system_information(RESULT _${key} QUERY ${key})
  cmake_print_variables(_${key})
endforeach()
message(STATUS "")
