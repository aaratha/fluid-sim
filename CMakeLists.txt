cmake_minimum_required(VERSION 3.14)
project(fluid-sim)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)  # Change to Debug if needed
endif()

# raylib
find_package(raylib QUIET)
if (NOT raylib_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        raylib
        GIT_REPOSITORY https://github.com/raysan5/raylib.git
        GIT_TAG 5.0
        GIT_SHALLOW 1
    )
    FetchContent_MakeAvailable(raylib)
endif()

# raylib-cpp
find_package(raylib_cpp QUIET)
if (NOT raylib_cpp_FOUND)
    if (NOT DEFINED RAYLIB_CPP_VERSION)
        set(RAYLIB_CPP_VERSION v5.0.2)
    endif()
    include(FetchContent)
    FetchContent_Declare(
        raylib_cpp
        GIT_REPOSITORY https://github.com/RobLoach/raylib-cpp.git
        GIT_TAG ${RAYLIB_CPP_VERSION}
    )
    FetchContent_MakeAvailable(raylib_cpp)
endif()

# raygui
include(FetchContent)
FetchContent_Declare(
    raygui
    GIT_REPOSITORY https://github.com/raysan5/raygui.git
    GIT_TAG 4.0
    GIT_SHALLOW 1
)
FetchContent_GetProperties(raygui)
if (NOT raygui_POPULATED)
    FetchContent_MakeAvailable(raygui)
endif()

# Include directories for the project
include_directories(
    ${CMAKE_SOURCE_DIR}/include  # Include your headers
    ${raylib_SOURCE_DIR}/src
    ${raygui_SOURCE_DIR}/src
)

# Source files
set(SOURCES
    src/main.cpp
    src/utils.cpp
    src/physics.cpp
)

# Main executable target
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC raylib raylib_cpp)


# Copy resources folder to the build directory
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/resources $<TARGET_FILE_DIR:${PROJECT_NAME}>/resources)


if(APPLE)
    target_link_libraries(${PROJECT_NAME} PUBLIC "-framework Cocoa" "-framework IOKit" "-framework CoreVideo")
endif()


# Finished setup!
