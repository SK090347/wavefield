cmake_minimum_required(VERSION 3.16)
project(wavefield_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(wavefield_cpp
  src/main.cpp
  src/solver.cpp
)
target_include_directories(wavefield_cpp PRIVATE src)

if(NOT MSVC)
  target_compile_options(wavefield_cpp PRIVATE -O3 -Wall -Wextra)
endif()

install(TARGETS wavefield_cpp RUNTIME DESTINATION bin)
