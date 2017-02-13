################################################################################
project( "TCO-CMAP2" C CXX )

cmake_minimum_required( VERSION 2.8.3 )

set( CMAKE_ECLIPSE_VERSION 3.8 )
set( _ECLIPSE_VERSION 3.8 )

set( CMAKE_CXX_COMPILER_ARG1 -std=c++11 ) ### for Eclipse's discovery extra arguments


#add_definitions( -O0 -ggdb -std=c++11 -Wall -pedantic -msse4.2 -march=westmere )
#add_definitions( -O2 -ggdb -std=c++11 -Wall -pedantic -mmmx -msse -msse2 -msse3 -s -pipe -Wno-sign-compare -march=native )
#add_definitions( -O2 -g -unroll-loops -std=c++11 -Wall -pedantic -msse4.2 -s -pipe -Wno-sign-compare -march=native )
add_definitions( -O3 -DNDEBUG -g -unroll-loops -std=c++11 -Wall -pedantic -s -pipe -Wno-sign-compare -msse4.2 -march=westmere )

#add_definitions( -O2 -ffast-math -fno-finite-math-only  -msse2 -march=native -std=c++11 -Wall -g -pedantic )
#add_definitions( -O2 -ftree-vectorize -ffast-math -fno-finite-math-only  -msse2 -march=native -std=c++11 -Wall -g -pedantic )
#add_definitions( -O3 -ffast-math -msse2 -march=native -std=c++11 -Wall -g -pedantic )

################################################################################

find_package( Boost REQUIRED COMPONENTS program_options )

################################################################################

include_directories(
    src
    cpplinq/CppLinq
    ${Boost_INCLUDE_DIR}
)

################################################################################

add_executable( main
    src/main.cpp
    src/parse_options.cpp
    src/ground_truth.cpp
    src/query_parser.cpp
    src/CMAPLib.cpp
)

target_link_libraries( main ${Boost_LIBRARIES} )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake-modules )

find_package( GTest REQUIRED )
find_package( GMock REQUIRED )
find_package( Threads REQUIRED )

include_directories(
    src
    ${GMOCK_INCLUDE_DIR}
    ${GTEST_INCLUDE_DIR}
)

add_executable( check
    src/test.cpp
)


target_link_libraries( check
    ${GTEST_BOTH_LIBRARIES}
    ${GMOCK_BOTH_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

add_executable( run_seq
    src/run_seq.cpp
    src/query_parser.cpp
)
target_link_libraries( run_seq
    ${CMAKE_THREAD_LIBS_INIT}
)

add_executable( run_conc
    src/run_conc.cpp
    src/query_parser.cpp
)
target_link_libraries( run_conc
    ${CMAKE_THREAD_LIBS_INIT}
)

################################################################################
