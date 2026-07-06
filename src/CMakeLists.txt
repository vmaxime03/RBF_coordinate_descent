cmake_minimum_required(VERSION 3.12...3.26)

project(ultimaille-examples LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(main main.cpp)

target_link_libraries(main  ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(main nlohmann_json::nlohmann_json)


add_executable(smoother smoother.cpp)
target_link_libraries(smoother ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(smoother nlohmann_json::nlohmann_json)

IF (NOT WIN32)
	target_link_libraries(main m)
	target_link_libraries(smoother m)
ENDIF()

set(OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../run_output")

file(MAKE_DIRECTORY "${OUTPUT_DIR}")
add_definitions(-DASSETS_INPUT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/assets/\")
add_definitions(-DOUTPUT_DIR=\"${OUTPUT_DIR}/\")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/script DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../assets")
add_definitions(-DINPUT_DIR=\"${INPUT_DIR}/\")

