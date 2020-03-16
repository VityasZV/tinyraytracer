cmake_minimum_required (VERSION 2.8)
project (raytracer)

include(CheckCXXCompilerFlag)

set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-pedantic")
enable_cxx_compiler_flag_if_supported("-std=c++17")
enable_cxx_compiler_flag_if_supported("-pthread")
enable_cxx_compiler_flag_if_supported("-O3")
#enable_cxx_compiler_flag_if_supported("-Ofast")



file(GLOB SOURCES *.h *.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} raytracing/raytracing.cpp picture/picture.cpp )



# STB
set(STB_LIB_DIR "${LIB_DIR}/stb_lib")
target_include_directories(${PROJECT_NAME} PRIVATE ${STB_LIB_DIR})

