cmake_minimum_required(VERSION 3.12...3.26)

project(ultimaille-examples LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(main main.cpp)

target_link_libraries(main  ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)

target_link_libraries(main nlohmann_json::nlohmann_json)

IF (NOT WIN32)
	target_link_libraries(main m)
ENDIF()

set(OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/output")

file(MAKE_DIRECTORY "${OUTPUT_DIR}")
add_definitions(-DASSETS_INPUT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/assets/\")
add_definitions(-DOUTPUT_DIR=\"${OUTPUT_DIR}/\")

# Copy asset files to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/script DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
