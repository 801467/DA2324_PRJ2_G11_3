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
