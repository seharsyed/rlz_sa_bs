cmake_minimum_required(VERSION 3.15)
project(sa_16mer_index)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3 -Wall -Wextra)

add_executable(clean_reference clean_reference.cpp)
add_executable(clean_inputs clean_inputs.cpp)

add_executable(pt16_rlz pt16_main.cpp)
target_include_directories(
    pt16_rlz
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../parser
)

add_executable(compare_rlz_cache compare_rlz_cache.cpp)
target_include_directories(
    compare_rlz_cache
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../parser
)

add_executable(pt16_build_time pt16_build_time.cpp)
target_include_directories(pt16_build_time PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../parser)