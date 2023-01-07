cmake_minimum_required(VERSION 3.21)

project(HEWES LANGUAGES C CXX VERSION 0.2.4)

# Print helper
include(CMakePrintHelpers)

# Print project name and version
cmake_print_variables(CMAKE_PROJECT_NAME)
cmake_print_variables(CMAKE_PROJECT_VERSION)

# Require C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check C++ compiler
cmake_print_variables(CMAKE_CXX_COMPILER_ID)
cmake_print_variables(CMAKE_CXX_COMPILER)
cmake_print_variables(CMAKE_CXX_STANDARD)

# If installed, do in build directory instead of system default
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/build")
cmake_print_variables(CMAKE_INSTALL_PREFIX)

# Warning on exotic system
if(NOT UNIX AND NOT WIN32)
    message(WARNING "Installation has only been tested on UNIX-like OS's \
    and Windows.")
endif()

# Built type
if(WIN32)  # Sundials only installs in Debug mode on Windows for whatever reason
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
    message(WARNING "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}, due to a \
    Sundials issue on Windows.")
else()  # Release is fastest
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Type of build (for SUNDIALS). Options are: \
    Debug, Release, RelWithDebInfo, MinSizeRel." FORCE)
    message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

# Dynamic or static libraries
if(WIN32)  # Otherwise problems with Sundials
    set(BUILD_SHARED_LIBS OFF)
else()
    option(BUILD_SHARED_LIBS "Build shared libraries." ON)
endif()

# Use of MPI
option(USE_MPI "Use of MPI" ON)
if(WIN32)  # MSPMI complies only with MPI 2.0 (adaptations in code)
    find_package(MPI 2.0 COMPONENTS CXX)
else()
    find_package(MPI 3.1 COMPONENTS CXX)
endif()
if(${MPI_FOUND} AND ${USE_MPI})
    message(STATUS "Using MPI.")
    cmake_print_variables(MPI_INCLUDE_PATH)
    cmake_print_variables(MPI_CXX_LIBRARIES)
    if(WIN32)  # Explicity point to MPI on Windows
        include_directories(SYSTEM ${MPI_INCLUDE_PATH})
        link_libraries(${MPI_CXX_LIBRARIES})
    endif()
elseif(NOT USE_MPI)
    message(STATUS "Not using MPI.")
endif()

# Use of OpenMP
option(USE_OPENMP "Use of OpenMP" ON)
find_package(OpenMP 4.5 COMPONENTS C CXX)  # C required for Sundials
if(${OpenMP_FOUND} AND ${USE_OPENMP})
    message(STATUS "Using OpenMP.")
    cmake_print_variables(OpenMP_C_INCLUDE_DIRS)
    cmake_print_variables(OpenMP_C_LIBRARIES)
    cmake_print_variables(OpenMP_CXX_INCLUDE_DIRS)
    cmake_print_variables(OpenMP_CXX_LIBRARIES)
elseif(NOT USE_OPENMP)
    message(STATUS "Not using OpenMP.")
endif()

# Find SUNDIALS
# Set path via `cmake -D SUNDIALS_DIR="/path/to/sundials/installation"` or here
if(NOT DEFINED SUNDIALS_DIR)
    set(SUNDIALS_DIR "~/Documents/sundials-6.3.0/INSTDIR_Rel")
endif()
# Search only in given path to avoid unintended package manager installation
find_package(SUNDIALS 6.0 NO_DEFAULT_PATH PATHS ${SUNDIALS_DIR})
if(${SUNDIALS_FOUND})
    message(STATUS "Found SUNDIALS in ${SUNDIALS_DIR}.")
else()
    message(WARNING "SUNDIALS not found. Fetching CVODE ...")
    include(FetchContent)
    if(${BUILD_SHARED_LIBS})
        set(BUILD_STATIC_LIBS OFF CACHE INTERNAL "")
        set(BUILD_SHARED_LIBS ON CACHE INTERNAL "")
    else()
        set(BUILD_STATIC_LIBS ON CACHE INTERNAL "")
        set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "")
    endif()
    if(${MPI_FOUND} AND ${USE_MPI})
        set(ENABLE_MPI ON CACHE INTERNAL "")
    endif()
    if(${OpenMP_FOUND} AND ${USE_OPENMP})
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
            GIT_TAG             v6.5.0)
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

# Compile options
if(MSVC)  # Otherwise 02 interference with RTC1
    set(CMAKE_CXX_FLAGS_DEBUG "")
endif()
target_compile_options(hewes PUBLIC -O2)
# Target host CPU
option(MARCH_NATIVE "Build with -march=native" OFF)
if(MARCH_NATIVE)
    add_compile_options(-march=native)
    message(STATUS "Building with -march=native.")
endif()

# Link Sundials, MPI, and OpenMP and the corresponding components
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
    HOSTNAME
    OS_NAME
    OS_RELEASE
    OS_VERSION
    OS_PLATFORM
    PROCESSOR_NAME
    PROCESSOR_DESCRIPTION
    NUMBER_OF_LOGICAL_CORES
    NUMBER_OF_PHYSICAL_CORES
    IS_64BIT
    TOTAL_VIRTUAL_MEMORY
    AVAILABLE_VIRTUAL_MEMORY
    TOTAL_PHYSICAL_MEMORY
    AVAILABLE_PHYSICAL_MEMORY
  )
  cmake_host_system_information(RESULT _${key} QUERY ${key})
  cmake_print_variables(_${key})
endforeach()
message(STATUS "")
