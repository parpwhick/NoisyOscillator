#cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
set (CMAKE_CXX_STANDARD 11y)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fuse-ld=gold -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g2 -Wall -Wextra -fno-gnu-unique -std=c++11 -fopenmp")
endif()


set (projectName NoisyOscillator)


	
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
	
set (Sources
	oscillator.cpp  
	Simulation.cpp  
	SimulationParameters.cpp
	Constants.h
	SimulationParameters.h
	Utilities.h
)

add_executable (${projectName} ${Sources})

target_link_libraries(${projectName} m)

