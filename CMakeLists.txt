cmake_minimum_required(VERSION 3.10)

project(catmodel)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add common library search paths for macOS (Homebrew/MacPorts)
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH 
        /usr/local 
        /opt/local 
        /opt/homebrew
    )
endif()

# Find required packages
find_package(PkgConfig REQUIRED)

# Try to find GiNaC using pkg-config first
pkg_check_modules(GINAC ginac)
if(NOT GINAC_FOUND)
    # Fallback: try to find manually
    find_path(GINAC_INCLUDE_DIR ginac/ginac.h)
    find_library(GINAC_LIBRARY ginac)
    if(GINAC_INCLUDE_DIR AND GINAC_LIBRARY)
        set(GINAC_FOUND TRUE)
        set(GINAC_INCLUDE_DIRS ${GINAC_INCLUDE_DIR})
        set(GINAC_LIBRARIES ${GINAC_LIBRARY})
    endif()
endif()

# Try to find CLN using pkg-config first
pkg_check_modules(CLN cln)
if(NOT CLN_FOUND)
    # Fallback: try to find manually
    find_path(CLN_INCLUDE_DIR cln/cln.h)
    find_library(CLN_LIBRARY cln)
    if(CLN_INCLUDE_DIR AND CLN_LIBRARY)
        set(CLN_FOUND TRUE)
        set(CLN_INCLUDE_DIRS ${CLN_INCLUDE_DIR})
        set(CLN_LIBRARIES ${CLN_LIBRARY})
    endif()
endif()

# Find GMP and GMPXX
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY gmp)
find_path(GMPXX_INCLUDE_DIR gmpxx.h)
find_library(GMPXX_LIBRARY gmpxx)

# Check if all libraries were found
if(NOT GINAC_FOUND)
    message(FATAL_ERROR "GiNaC library not found")
endif()

if(NOT CLN_FOUND)
    message(FATAL_ERROR "CLN library not found")
endif()

if(NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "GMP or GMPXX library not found")
endif()

# Create executable
add_executable(catmodel catmodel.cc)

# Link libraries (note: order matters, GiNaC depends on CLN and GMP)
target_link_libraries(catmodel 
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
    ${GMPXX_LIBRARY}
    ${GMP_LIBRARY}
)

# Include directories
target_include_directories(catmodel PRIVATE
    ${GINAC_INCLUDE_DIRS}
    ${CLN_INCLUDE_DIRS}
    ${GMP_INCLUDE_DIR}
    ${GMPXX_INCLUDE_DIR}
)

# Compiler flags
target_compile_options(catmodel PRIVATE
    ${GINAC_CFLAGS_OTHER}
    ${CLN_CFLAGS_OTHER}
)