cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(VesselExtraction)

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

SET(Sources src/typedefs.h src/bloodvesselpattern.h src/bloodvesselpattern.cpp src/segmentation.cpp)

ADD_EXECUTABLE(VesselExtraction ${Sources} )

TARGET_LINK_LIBRARIES(VesselExtraction ITKCommon ITKIO)
SET_TARGET_PROPERTIES(VesselExtraction PROPERTIES LINK_FLAGS "-openmp")