CMAKE_MINIMUM_REQUIRED(VERSION 3.13.5)

find_package(OpenMP REQUIRED)

set(sfa_SOURCES 
    src/sfa.cpp
    src/utils.cpp)

add_library(sfa SHARED ${sfa_SOURCES})
target_include_directories(sfa PUBLIC ${CMAKE_CURRENT_BINARY_DIR} include)
target_compile_features(sfa PRIVATE cxx_std_17)
target_link_libraries(sfa PUBLIC OpenMP::OpenMP_CXX)
set(CMAKE_CXX_FLAGS "-O3 -lm -Wall -Wextra -ffast-math")




