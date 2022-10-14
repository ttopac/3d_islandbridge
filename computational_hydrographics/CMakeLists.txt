cmake_minimum_required(VERSION 2.6)
project(ig)

add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_ZERO)

#SET(LIBIGL_USE_STATIC_LIBRARY 1)
#add_definitions(-DIGL_STATIC_LIBRARY)

message(STATUS "CMAKE_CXX_FLAGS: ${PROJECT_SOURCE_DIR}/../libigl/tutorial/cmake")

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../libigl/tutorial/cmake)

IF(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") ### Enable parallel compilation for Visual Studio
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
ELSE(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../" )
ENDIF(MSVC)

if(APPLE)
	set(CMAKE_CXX_LINK_FLAGS "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
endif (APPLE) #APPLE



find_package(EMBREE REQUIRED)
include_directories(${EMBREE_INCLUDE_DIRS})

find_package(EIGEN REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})

#Compile libigl (uncomment to use static library)
#add_definitions(-DIGL_NO_OPENGL)
#add_definitions(-DIGL_STATIC_LIBRARY)
#add_subdirectory("../libigl/optional" "libigl")
include_directories("../libigl/include")

#Compile nanogui
set(NANOGUI_BUILD_PYTHON OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_SHARED OFF CACHE BOOL " " FORCE)
add_subdirectory("../libigl/external/nanogui/" "nanogui")
include_directories("../libigl/external/nanogui/include")
include_directories("../libigl/external/nanogui/ext/nanovg/src")
include_directories("../libigl/external/nanogui/ext/glfw/include")
list(APPEND SHARED_LIBRARIES "nanogui" "glfw")

#Compile libpng
option(PNG_SHARED "Build shared lib" OFF)
add_subdirectory("../libigl/external/libpng/" "libpng")
include_directories("../libigl/external/libpng/")

message(STATUS "Linking with: ${LIBIGL_LIBRARIES}")
list(APPEND SHARED_LIBRARIES ${LIBIGL_LIBRARIES})

# This is necessary to work around the static linking order dependencies on linux
if(UNIX AND NOT APPLE)
  set(SHARED_LIBRARIES ${SHARED_LIBRARIES} ${SHARED_LIBRARIES})
endif(UNIX AND NOT APPLE)

# add_definitions(-DIGL_STATIC_LIBRARY)

file(GLOB SOURCES_SIM
    "${PROJECT_SOURCE_DIR}/src/simulation/*.h"
    "${PROJECT_SOURCE_DIR}/src/simulation/*.cpp"
)

LIST(APPEND SOURCES_SIM ${YIMG_SOURCES})

message(STATUS ${SHARED_LIBRARIES})

add_executable(${PROJECT_NAME}_sim ${SOURCES_SIM} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_sim ${SHARED_LIBRARIES} png16_static ${EMBREE_LIBRARIES})
