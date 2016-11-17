cmake_minimum_required (VERSION 2.8)
project (AGO)

set( CMAKE_VERBOSE_MAKEFILE ON )
set( OPTIMIZE -O3 )
set( DEBUG "-g")
set( WARN "-Wno-ignored-attributes -Wno-multichar -Wreturn-type -Wno-terminate" )

set( CMAKE_CXX_FLAGS "-std=c++17 -ftemplate-backtrace-limit=0 ${OPTIMIZE} ${DEBUG} ${WARN} " )

set( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/local/lib/" )

if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    execute_process( COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION )
    if( NOT (GCC_VERSION VERSION_GREATER 5.0 OR GCC_VERSION VERSION_EQUAL 5.0) )
        message( FATAL_ERROR "${PROJECT_NAME} requires g++ 5.0 OR greater." )
    endif()
else()
    message( FATAL_ERROR "Your C++ compiler does not support C++14." )
endif()

set( Boost_USE_STATIC_LIBS    OFF )
set( Boost_USE_MULTITHREADED  ON  )
set( Boost_USE_STATIC_RUNTIME OFF )
find_package( Boost COMPONENTS 
    date_time 
    program_options 
    filesystem system 
    serialization 
    regex 
    thread 
    iostreams )
find_package( Threads )
find_package( ZLIB )
find_package( CURL )
find_package( GTest REQUIRED )

if( Boost_FOUND AND Threads_FOUND AND ZLIB_FOUND )

    include_directories( ${CMAKE_SOURCE_DIR}/include )
    include_directories( ${GTEST_INCLUDE_DIRS} )

    link_directories( ${CMAKE_SOURCE_DIR}/lib ${CMAKE_SOURCE_DIR}/lib/* )
    set( EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin )
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/ )
    execute_process( COMMAND ln -sf ../bin )
    execute_process( COMMAND ln -sf ../src )
    execute_process( COMMAND ln -sf ../include )

    add_executable( pipeline_builder src/application/pipeline_builder/main.cpp )
    target_link_libraries( pipeline_builder ${Boost_LIBRARIES} ${CURL_LIBRARIES} ) 

endif()
