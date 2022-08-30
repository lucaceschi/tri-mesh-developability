cmake_minimum_required(VERSION 3.20.5)

project(developability)

set(CMAKE_CXX_STANDARD 14)

set(VCGLIB_DIR vcglib)
set(SRC_DIR src)
set(INCL_DIR include)
set(SOURCES
	${SRC_DIR}/main.cpp
	${SRC_DIR}/mesh.cpp
	${SRC_DIR}/energy.cpp
	${SRC_DIR}/energy_grad.cpp
	${SRC_DIR}/opt.cpp
	${SRC_DIR}/remeshing.cpp
	${INCL_DIR}/mesh.hpp
	${INCL_DIR}/energy.hpp
	${INCL_DIR}/energy_grad.hpp
	${INCL_DIR}/opt.hpp
	${INCL_DIR}/remeshing.hpp
	${VCGLIB_DIR}/wrap/ply/plylib.cpp
)

add_subdirectory(${VCGLIB_DIR})

add_executable(develop ${SOURCES})
target_include_directories(develop PRIVATE ${INCL_DIR})
target_link_libraries(develop
	PUBLIC vcglib
)