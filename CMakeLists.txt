cmake_minimum_required(VERSION 3.20.5)

project(developability)

set(CMAKE_CXX_STANDARD 14)

set(VCGLIB_DIR vcglib)
set(SRC_DIR src)
set(INCL_DIR include)
set(SOURCES
	${SRC_DIR}/main.cpp
	${INCL_DIR}/mesh.hpp
	${VCGLIB_DIR}/wrap/ply/plylib.cpp
)

add_subdirectory(${VCGLIB_DIR})

add_executable(out ${SOURCES})
target_include_directories(out PRIVATE ${INCL_DIR})
target_link_libraries(out
	PUBLIC vcglib
)
