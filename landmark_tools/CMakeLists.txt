cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(Converter)

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
INCLUDE(${CMAKE_ROOT}/Modules/FindVTK.cmake)
IF(USE_VTK_FILE)
  INCLUDE(${USE_VTK_FILE})
ELSE(USE_VTK_FILE)
  MESSAGE(FATAL_ERROR
          "Cannot build without VTK.  Please set VTK_DIR.")
ENDIF(USE_VTK_FILE)

#Find GSL
FIND_PACKAGE(GSL REQUIRED)

#CMake Optimizing
SET(CMAKE_BUILD_TYPE RELEASE )
SET(CMAKE_CXX_FLAGS_DEBUG "-g -wd1268")
SET(CMAKE_C_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -openmp")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -openmp")
SET(CMAKE_VERBOSE_MAKEFILE "TRUE")

#Sources
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

SET(Sources1 converter.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp)
SET(Sources2 landmark_correspondence.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp)
SET(Sources3 tripleLabelingAnalysis.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp)
SET(Sources4 joinCommonLandmarks.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp)

ADD_EXECUTABLE(Converter ${Sources1} )
ADD_EXECUTABLE(LandmarkCorrespondence ${Sources2} )
ADD_EXECUTABLE(TripleLabelingCorrespondence ${Sources3} )
ADD_EXECUTABLE(JoinCommonLandmarks ${Sources4} )

TARGET_LINK_LIBRARIES(Converter ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(Converter PROPERTIES LINK_FLAGS "-openmp")

TARGET_LINK_LIBRARIES(LandmarkCorrespondence ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(LandmarkCorrespondence PROPERTIES LINK_FLAGS "-openmp")

TARGET_LINK_LIBRARIES(TripleLabelingCorrespondence ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(TripleLabelingCorrespondence PROPERTIES LINK_FLAGS "-openmp")

TARGET_LINK_LIBRARIES(JoinCommonLandmarks ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(JoinCommonLandmarks PROPERTIES LINK_FLAGS "-openmp")