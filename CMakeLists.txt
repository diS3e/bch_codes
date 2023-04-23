cmake_minimum_required(VERSION 3.24)
project(bch_codes)

set(CMAKE_CXX_STANDARD 17)

find_package(GTest REQUIRED)

target_link_libraries(tests GTest::gtest GTest::gtest_main)

add_executable(bch_codes main.cpp tests.cpp)
