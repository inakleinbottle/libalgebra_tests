cmake_minimum_required(VERSION 3.10)

option(USE_FLATMAP "Use Boost flatmap instead of std::map" OFF)


project(LibAlgebraUnitTests VERSION 1.0)

find_package(UnitTest++ REQUIRED)

find_package(Boost COMPONENTS filesystem iostreams REQUIRED)
if (USE_FLATMAP)
    find_package(Boost COMPONENTS container filesystem iostreams REQUIRED)
endif()



add_compile_definitions(DENSE_TRANSITION=2<<5)

find_package(benchmark REQUIRED)


include_directories(libalgebra)

file(GLOB SOURCES "LibAlgebraUnitTests/*.cpp")

configure_file(libalgebra_unittests.h.in unittest_config.h)


#add_executable(test "LibAlgebraUnitTests/LibAlgebraUnitTests.cpp" "LibAlgebraUnitTests/AlgebraFunctionsTests.cpp" "LibAlgebraUnitTests/memfile.cpp" "LibAlgebraUnitTests/makebm.cpp")
add_executable(test ${SOURCES})
target_include_directories(test PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(test UnitTest++ gmp ${Boost_LIBRARIES} pthread)
set_property(TARGET test PROPERTY CXX_STANDARD 11)


file(GLOB BENCH_SOURCES "benchmarks/*.cpp")

add_executable(benchmarks ${BENCH_SOURCES})
target_include_directories(benchmarks PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(benchmarks benchmark::benchmark UnitTest++ gmp ${Boost_LIBRARIES} pthread)
set_property(TARGET test PROPERTY CXX_STANDARD 11)



file(GLOB SVD_SOURCES "svd/*.cpp")

add_executable(time_vec_iteration "${SVD_SOURCES}")
target_include_directories(time_vec_iteration PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(time_vec_iteration benchmark::benchmark gmp ${Boost_LIBRARIES} pthread)
set_property(TARGET test PROPERTY CXX_STANDARD 11)