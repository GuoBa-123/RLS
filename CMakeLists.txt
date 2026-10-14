cmake_minimum_required(VERSION 3.20)
project(rls VERSION 0.1.0 LANGUAGES CXX)

# header-only 库
add_library(rls INTERFACE)
add_library(rls::rls ALIAS rls)

target_include_directories(rls INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(rls INTERFACE cxx_std_17)

