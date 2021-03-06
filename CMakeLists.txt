cmake_minimum_required (VERSION 2.8)
project (AGO)

set( CMAKE_VERBOSE_MAKEFILE ON )
set( LOGGER_SETTING -DLOGGER_OPEN_CHANNEL="0,1,2,4")
set( OPTIMIZE -Ofast )
# set( DEBUG "-g" )
set( GUI_MODEL -DCPT_MONITOR_TGUI )
set( WARN "-Wno-multichar -Wno-return-type -Wno-write-strings -Wno-format" )
set( SINGLE_CPP "-DSINGLE_CPP" )

set( CMAKE_CXX_FLAGS "-std=c++17 -ftemplate-backtrace-limit=0 ${LOGGER_SETTING} ${OPTIMIZE} ${DEBUG} ${SINGLE_CPP} ${WARN} ${GUI_MODEL} " )

set( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/local/lib/" )

if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    execute_process( COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION )
    if( NOT (GCC_VERSION VERSION_GREATER 5.0 OR GCC_VERSION VERSION_EQUAL 5.0) )
        message( FATAL_ERROR "${PROJECT_NAME} requires g++ 5.0 OR greater." )
    endif()
else()
    message( FATAL_ERROR "Your C++ compiler does not support C++14." )
endif()

include( ${CMAKE_ROOT}/Modules/FindBoost.cmake )
include( ${CMAKE_ROOT}/Modules/FindGTest.cmake )

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
find_package( ZLIB    )
find_package( CURL    )
find_package( GTest   )

if( Boost_FOUND AND Threads_FOUND AND ZLIB_FOUND )

    include_directories( ${CMAKE_SOURCE_DIR}/include )
    include_directories( ${CMAKE_SOURCE_DIR}/lib/CCD/include )
    include_directories( ${CMAKE_SOURCE_DIR}/lib/CPT/include )
    include_directories( ${CMAKE_SOURCE_DIR}/lib/slam_quant_seq/include )
    include_directories( ${CMAKE_SOURCE_DIR}/lib/Tailor/include )
    include_directories( ${GTEST_INCLUDE_DIRS} )

    link_directories( ${CMAKE_SOURCE_DIR}/lib ${CMAKE_SOURCE_DIR}/lib/* )

    set( EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin )
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/ )

    add_executable( pipeline_builder src/application/pipeline_builder/main.cpp )
        target_link_libraries( pipeline_builder
        ${Boost_LIBRARIES} ${CURL_LIBRARIES}
        pthread thread_pool ) 

    # add_executable( md_sam_test unit_test/format/md_sam.cpp )
    #     target_link_libraries( md_sam_test
    #     ${Boost_LIBRARIES} ${CURL_LIBRARIES}
    #     gtest gtest_main pthread thread_pool boost_program_options ) 

    # add_executable( tailor_test unit_test/aligner/tailor.cpp )
    #     target_link_libraries( tailor_test
    #     ${Boost_LIBRARIES} ${CURL_LIBRARIES}
    #     gtest gtest_main pthread thread_pool boost_program_options ) 

    # add_executable( annotator_test unit_test/annotator/annotator.cpp )
    #     target_link_libraries( annotator_test
    #     ${Boost_LIBRARIES} ${CURL_LIBRARIES}
    #     gtest gtest_main pthread thread_pool ) 

    # add_executable( rnafold src/application/rnafold/main.cpp )
    #     target_link_libraries( rnafold ${Boost_LIBRARIES} )

    # add_executable( gtf2bed src/application/gtf2bed/main.cpp )
    #     target_link_libraries( gtf2bed ${Boost_LIBRARIES} )

    # add_executable( arcReader src/application/arc_reader/main.cpp )
    #     target_link_libraries( arcReader ${Boost_LIBRARIES} )

    # add_executable( getCounts src/application/get_counts/main.cpp )
    #     target_link_libraries( getCounts ${Boost_LIBRARIES} )

    # add_executable( makeRank src/application/make_rank/main.cpp )
    #     target_link_libraries( makeRank ${Boost_LIBRARIES} )

endif()
