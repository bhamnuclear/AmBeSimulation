#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(ambe)

#----------------------------------------------------------------------------
# check MPI package...
find_package(MPI REQUIRED)

#----------------------------------------------------------------------------
# Find ROOT package, required by TGraphs and ROOT specific objects
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIRS})
message(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})
#link_directories(${ROOT_LIBRARIES})
message(${ROOT_LIBRARY_DIR})

# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

find_package(G4mpi REQUIRED)

#----------------------------------------------------------------------------
# Setup ROOT include directories and compile definitions
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
link_directories(${ROOT_LIBRARY_DIR})
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${G4mpi_INCLUDE_DIR})



# We presume the existence of three directories
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB_RECURSE macros RELATIVE ${PROJECT_SOURCE_DIR} macros/*.mac)
file(GLOB_RECURSE xSections RELATIVE ${PROJECT_SOURCE_DIR} AmBeData/*)

# Enable macros for out-of-source build
foreach(_file ${macros})
  configure_file(
    ${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
    )
endforeach()

foreach(_file ${xSections})
  configure_file(
    ${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
    )
endforeach()

#Copy script
configure_file(
    ${PROJECT_SOURCE_DIR}/MakePlots.cc
    ${PROJECT_BINARY_DIR}/MakePlots.cc
    COPYONLY
    )



add_executable(ambe main.cc ${sources} ${headers})
target_link_libraries(ambe ${G4mpi_LIBRARIES} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})


install(TARGETS ambe DESTINATION bin)
