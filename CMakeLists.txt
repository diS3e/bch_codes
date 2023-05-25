cmake_minimum_required(VERSION 3.24)
project(bch_codes)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(lib)
include_directories(lib/googletest/include)
add_executable(tests tests.cpp bch_codes.h bch_codes.cpp galois_field.h galois_field.cpp chase_decoder.h random.h)
add_executable(main main.cpp bch_codes.h bch_codes.cpp galois_field.h galois_field.cpp chase_decoder.h random.h)
add_executable(experiment experiment.cpp bch_codes.h bch_codes.cpp galois_field.h galois_field.cpp chase_decoder.h random.h)
target_link_libraries(tests GTest::gtest GTest::gtest_main)
