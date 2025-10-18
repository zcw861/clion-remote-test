cmake_minimum_required(VERSION 3.29)
project(test)

set(CMAKE_CXX_STANDARD 20)

# set(test src/main.cpp)
# aux_source_directory(${PROJECT_SOURCE_DIR} src/test)
file(GLOB test ${CMAKE_CURRENT_SOURCE_DIR}/CMake_test/src/*.cpp)

set(SOURCES ${test} ${main})

include_directories(include)

#add_library(calc SHARED ${test})
add_library(calc STATIC ${test})

add_executable(app main.cpp)

target_link_libraries(app calc)