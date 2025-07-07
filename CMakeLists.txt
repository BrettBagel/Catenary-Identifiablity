cmake_minimum_required(VERSION 3.16)
project(ginac_project)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# SET(CMAKE_CXX_COMPILER "/usr/local/bin/g++-14" CACHE STRING "C++ compiler" FORCE)

# Find required packages
find_package(PkgConfig REQUIRED)

# Find CLN library
pkg_check_modules(CLN REQUIRED cln)

# Find GiNaC library
pkg_check_modules(GINAC REQUIRED ginac)

# pkg_check_modules(GMP )

# Add executable
add_executable(${PROJECT_NAME} catmodel.cc)

# Link libraries
target_link_libraries(${PROJECT_NAME} 
    ${CLN_LIBRARIES}
    ${GINAC_LIBRARIES}
)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CLN_INCLUDE_DIRS}
    ${GINAC_INCLUDE_DIRS}
)

# Compile definitions
target_compile_definitions(${PROJECT_NAME} PRIVATE
    ${CLN_CFLAGS_OTHER}
    ${GINAC_CFLAGS_OTHER}
)

# Set library directories
target_link_directories(${PROJECT_NAME} PRIVATE
    ${CLN_LIBRARY_DIRS}
    ${GINAC_LIBRARY_DIRS}
)