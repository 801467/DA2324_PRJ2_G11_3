cmake_minimum_required(VERSION 3.27)
project(DA2324_PRJ2_G11_3)

set(CMAKE_CXX_STANDARD 14)

include_directories(classes)
include_directories(functionalities)
include_directories(terminal)
include_directories(utils)

add_executable(DA2324_PRJ2_G11_3
        classes/Graph.h
        classes/TSPGraph.cpp
        functionalities/Functionality.cpp
        terminal/TerminalFlow.cpp
        utils/FileReader.cpp
        main.cpp
        utils/haversine.cpp
)

# Doxygen Build
find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(BUILD_DOC_DIR "${CMAKE_SOURCE_DIR}/docs/output")
    if(NOT EXISTS "${BUILD_DOC_DIR}")
        file(MAKE_DIRECTORY "${BUILD_DOC_DIR}")
    endif()#

    set(DOXYGEN_IN "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
    set(DOXYGEN_OUT "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
    configure_file("${DOXYGEN_IN}" "${DOXYGEN_OUT}" @ONLY)

    message("Doxygen build started")
    add_custom_target(Doxygen ALL
            COMMAND "${DOXYGEN_EXECUTABLE}" "${DOXYGEN_OUT}"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
else(DOXYGEN_FOUND)
    message("Doxygen needs to be installed to generate the documentation.")
endif(DOXYGEN_FOUND)

target_link_libraries(DA2324_PRJ2_G11_3)
