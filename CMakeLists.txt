cmake_minimum_required(VERSION 2.8.5)
set (CMAKE_CXX_STANDARD 14)
# set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fuse-ld=gold")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g2 -Wall -Wextra -fno-gnu-unique -std=c++11")
endif()


set (projectName NoisyOscillator)


include_directories(Eigen)	
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
	
set (Sources
	oscillator.cpp  
	Simulation.cpp  
	Simulation.hpp
	SimulationParameters.cpp
	Constants.h
	SimulationParameters.h
	Utilities.cpp
	Utilities.h
)

add_executable (${projectName} ${Sources})
target_link_libraries(${projectName} m)
