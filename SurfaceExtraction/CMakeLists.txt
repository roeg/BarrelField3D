cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(SurfaceExtraction)

#CMake Optimizing
#SET(CMAKE_BUILD_TYPE Release )
#SET(CMAKE_CXX_FLAGS_Release "-O2")
#SET(CMAKE_C_FLAGS_Release "-O2")

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Find VTK
#INCLUDE(${CMAKE_ROOT}/Modules/FindVTK.cmake)
#IF(USE_VTK_FILE)
#  INCLUDE(${USE_VTK_FILE})
#ELSE(USE_VTK_FILE)
#  MESSAGE(FATAL_ERROR
#          "Cannot build without VTK.  Please set VTK_DIR.")
#ENDIF(USE_VTK_FILE)

#LINK_LIBRARIES(ITKBasicFilters ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid)
LINK_LIBRARIES(ITKBasicFilters ITKIO ITKAlgorithms)

#CMake Optimizing
SET(CMAKE_BUILD_TYPE RELEASE )
SET(CMAKE_CXX_FLAGS_DEBUG "-g -wd1268")
SET(CMAKE_C_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -openmp")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -openmp")
SET(CMAKE_VERBOSE_MAKEFILE "TRUE")

#Sources
INCLUDE_DIRECTORIES(src)

#SET(Sources src/typedefs.h src/surface_extract.cpp src/segmentation.h src/segmentation.cpp src/basics.h src/basics.cpp src/bloodvesselpattern.h src/bloodvesselpattern.cpp src/amiraReader.h src/amiraReader.cpp src/spatialGraph.h src/spatialGraph.cpp src/section.h src/section.cpp)
SET(Sources src/typedefs.h src/surface_extract.cpp src/segmentation.h src/segmentation.cpp src/basics.h src/basics.cpp src/bloodvesselpattern.h src/bloodvesselpattern.cpp src/section.h src/section.cpp)

ADD_EXECUTABLE(SurfaceExtraction ${Sources} )

TARGET_LINK_LIBRARIES(SurfaceExtraction ITKCommon ITKIO)
SET_TARGET_PROPERTIES(SurfaceExtraction PROPERTIES LINK_FLAGS "-openmp")
