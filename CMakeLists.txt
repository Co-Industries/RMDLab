cmake_minimum_required(VERSION 3.5)
project(RMD)
include_directories(${CMAKE_BINARY_DIR})
find_package(Git)
add_custom_target(version
    ${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/Version.h.in
                        -D DST=${CMAKE_BINARY_DIR}/Version.h
                        -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                        -P ${CMAKE_SOURCE_DIR}/GenerateVersionHeader.cmake
)    
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules/" ${CMAKE_MODULE_PATH})

# Add MagnumIntegration as a subproject
add_subdirectory(magnum-integration EXCLUDE_FROM_ALL)

# Add Corrade as a subproject
add_subdirectory(corrade EXCLUDE_FROM_ALL)

# Add Magnum as a subproject, enable Sdl2Application
set(CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/SDL2-2.0.10 ${CMAKE_PREFIX_PATH})
set(MAGNUM_WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)
add_subdirectory(magnum EXCLUDE_FROM_ALL)

add_subdirectory(src)
