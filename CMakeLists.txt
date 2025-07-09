cmake_minimum_required(VERSION 3.10)

project(catmodel)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add platform-specific library search paths
if(APPLE)
    # macOS with Homebrew
    execute_process(COMMAND brew --prefix OUTPUT_VARIABLE HOMEBREW_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(HOMEBREW_PREFIX)
        list(APPEND CMAKE_PREFIX_PATH ${HOMEBREW_PREFIX})
    endif()
    # Also check common macOS paths
    list(APPEND CMAKE_PREFIX_PATH /usr/local /opt/local /opt/homebrew)
endif()

# Find required packages
find_package(PkgConfig REQUIRED)

# Function to find library with fallback
function(find_library_with_fallback LIB_NAME HEADER_NAME)
    string(TOUPPER ${LIB_NAME} LIB_UPPER)
    
    # Try pkg-config first
    pkg_check_modules(${LIB_UPPER} ${LIB_NAME})
    
    if(NOT ${LIB_UPPER}_FOUND)
        # Fallback to manual search
        find_path(${LIB_UPPER}_INCLUDE_DIR ${HEADER_NAME})
        find_library(${LIB_UPPER}_LIBRARY ${LIB_NAME})
        
        if(${LIB_UPPER}_INCLUDE_DIR AND ${LIB_UPPER}_LIBRARY)
            set(${LIB_UPPER}_FOUND TRUE PARENT_SCOPE)
            set(${LIB_UPPER}_INCLUDE_DIRS ${${LIB_UPPER}_INCLUDE_DIR} PARENT_SCOPE)
            set(${LIB_UPPER}_LIBRARIES ${${LIB_UPPER}_LIBRARY} PARENT_SCOPE)
            set(${LIB_UPPER}_CFLAGS_OTHER "" PARENT_SCOPE)
        endif()
    endif()
endfunction()

# Find libraries
find_library_with_fallback(ginac ginac/ginac.h)
find_library_with_fallback(cln cln/cln.h)

# Find GMP and GMPXX
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY gmp)
find_path(GMPXX_INCLUDE_DIR gmpxx.h)
find_library(GMPXX_LIBRARY gmpxx)

# Check if all libraries were found
if(NOT GINAC_FOUND)
    message(FATAL_ERROR "GiNaC library not found. Please install libginac-dev (Ubuntu/Debian) or ginac-devel (Fedora/RHEL) or ginac (Homebrew/MacPorts)")
endif()

if(NOT CLN_FOUND)
    message(FATAL_ERROR "CLN library not found. Please install libcln-dev (Ubuntu/Debian) or cln-devel (Fedora/RHEL) or cln (Homebrew/MacPorts)")
endif()

if(NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "GMP or GMPXX library not found. Please install libgmp-dev (Ubuntu/Debian) or gmp-devel (Fedora/RHEL) or gmp (Homebrew/MacPorts)")
endif()

# Create executable
add_executable(catmodel catmodel.cc)

# Link libraries (order matters: GiNaC depends on CLN and GMP)
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

# Print found libraries for debugging
message(STATUS "GiNaC: ${GINAC_LIBRARIES}")
message(STATUS "CLN: ${CLN_LIBRARIES}")
message(STATUS "GMP: ${GMP_LIBRARY}")
message(STATUS "GMPXX: ${GMPXX_LIBRARY}")