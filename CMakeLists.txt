cmake_minimum_required(VERSION 3.16)

project(Spectrum_Init)

set(CMAKE_CXX_STANDARD 20)

add_executable(Spectrum_Init src/main.cpp)
target_link_libraries(Spectrum_Init fftw3)