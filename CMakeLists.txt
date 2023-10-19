CMAKE_MINIMUM_REQUIRED(VERSION 3.13.5)

project(pysfa)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "-O3 -lm -Wall -Wextra -ffast-math -g -march=native")
set(pysfa_SOURCES 
    src/pysfa.cpp
    src/utils.cpp)

if (TARGET propag-plugins)
    list(APPEND pysfa_SOURCES src/harmonic_pysfa.cpp)
    add_library(pysfa ${pysfa_SOURCES})
    target_link_libraries(pysfa PUBLIC propag)
    add_dependencies(pysfa propag)
else()
    add_library(pysfa SHARED ${pysfa_SOURCES})
endif()
target_include_directories(pysfa PUBLIC ${CMAKE_CURRENT_BINARY_DIR} include)
target_compile_features(pysfa PRIVATE cxx_std_17)
target_link_libraries(pysfa PUBLIC OpenMP::OpenMP_CXX)

