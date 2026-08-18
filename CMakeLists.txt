cmake_minimum_required(VERSION 3.15)
project(telox_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make Conda-provided CMake packages (re2, Abseil and llvm-openmp) available
# during pip's isolated editable build without replacing its pybind11 prefix.
if(DEFINED ENV{CONDA_PREFIX} AND EXISTS "$ENV{CONDA_PREFIX}")
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(bloom_core src/mnode.cpp)
target_compile_options(bloom_core PRIVATE -O3)
install(TARGETS bloom_core DESTINATION telox)

pybind11_add_module(telohmm src/telohmm.cpp)
if(NOT MSVC) 
    target_compile_options(telohmm PRIVATE -O3 -ffast-math -fno-strict-aliasing)
endif()
install(TARGETS telohmm DESTINATION telox)

pybind11_add_module(motiflev src/motiflev.cpp)
if(NOT MSVC) 
    target_compile_options(motiflev PRIVATE -O3 -ffast-math -fno-strict-aliasing)
else()
    target_compile_options(motiflev PRIVATE /O2)
endif()
install(TARGETS motiflev DESTINATION telox)

find_package(ZLIB REQUIRED)
find_package(OpenMP)
find_package(re2 REQUIRED)

add_executable(telogrep src/telogrep.cpp)
target_link_libraries(telogrep PRIVATE ZLIB::ZLIB re2::re2)
if(OpenMP_CXX_FOUND)
    target_link_libraries(telogrep PRIVATE OpenMP::OpenMP_CXX)
endif()

# A pip editable install puts telogrep in <prefix>/bin while Conda libraries
# remain in <prefix>/lib.  Preserve that relationship at runtime on macOS.
if(APPLE)
    set_target_properties(telogrep PROPERTIES
        INSTALL_RPATH "@loader_path/../lib"
    )
endif()

if(SKBUILD)
    install(TARGETS telogrep DESTINATION "${SKBUILD_SCRIPTS_DIR}")
else()
    install(TARGETS telogrep DESTINATION bin)
endif()
