cmake_minimum_required(VERSION 3.10)

option(PROFILER "Enable profiling features" OFF)


project(LibAlgebraUnitTests VERSION 1.0)

find_package(UnitTest++ REQUIRED)

find_package(Boost COMPONENTS filesystem iostreams serialization REQUIRED)


if (PROFILER)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -pg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -pg")
endif()

add_compile_definitions(DENSE_TRANSITION=2<<5)

find_package(benchmark REQUIRED)


include_directories(libalgebra)

file(GLOB SOURCES    "LibAlgebraUnitTests/*.cpp")

configure_file(libalgebra_unittests.h.in unittest_config.h)


#add_executable(test "LibAlgebraUnitTests/LibAlgebraUnitTests.cpp" "LibAlgebraUnitTests/AlgebraFunctionsTests.cpp" "LibAlgebraUnitTests/memfile.cpp" "LibAlgebraUnitTests/makebm.cpp")
add_executable(test ${SOURCES})
target_include_directories(test PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(test UnitTest++ ${Boost_LIBRARIES} pthread gmpxx gmp)
set_property(TARGET test PROPERTY CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" -pg")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

file(GLOB BENCH_SOURCES "benchmarks/*.cpp")

add_executable(benchmarks ${BENCH_SOURCES} benchmarks/time_index_of_key.cpp.d)
target_include_directories(benchmarks PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(benchmarks benchmark::benchmark UnitTest++ ${Boost_LIBRARIES} pthread gmpxx gmp)
set_property(TARGET test PROPERTY CXX_STANDARD 11)



file(GLOB SVD_SOURCES "svd/*.cpp")

add_executable(time_vec_iteration "${SVD_SOURCES}")
target_include_directories(time_vec_iteration PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(time_vec_iteration benchmark::benchmark ${Boost_LIBRARIES} pthread gmpxx gmp )
#set_property(TARGET test PROPERTY CXX_STANDARD 11)