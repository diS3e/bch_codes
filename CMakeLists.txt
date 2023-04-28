cmake_minimum_required(VERSION 3.24)
project(bch_codes)

set(CMAKE_CXX_STANDARD 17)

find_package(GTest REQUIRED)


add_executable(tests tests.cpp bch_codes.h bch_codes.cpp)
target_link_libraries(tests GTest::gtest GTest::gtest_main)
