cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(AverageBarrelField)

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

#LINK_LIBRARIES(ITKBasicFilters ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid)

#CMake Optimizing
SET(CMAKE_BUILD_TYPE RELEASE )
SET(CMAKE_CXX_FLAGS_DEBUG "-g -wd1268")
SET(CMAKE_C_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -openmp")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -openmp")
SET(CMAKE_VERBOSE_MAKEFILE "TRUE")

#Sources
# INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(src ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

SET(Sources src/main.cpp src/registration.h src/registration.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp)
#SET(Sources2 src/barrel_contours.cpp src/voronoiRegions.h src/voronoiRegions.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp ../common/imageIO.h ../common/imageIO.cpp)
SET(Sources3 src/contour_centroid.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp)
SET(AxisFieldSrcs src/axis_vectorfield.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp)

ADD_EXECUTABLE(AverageBarrelField ${Sources} )
#ADD_EXECUTABLE(BarrelContours ${Sources2} )
# ADD_EXECUTABLE(LandmarkCentroids ${Sources3} )
ADD_EXECUTABLE(AxisVectorField ${AxisFieldSrcs} )

TARGET_LINK_LIBRARIES(AverageBarrelField ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(AverageBarrelField PROPERTIES LINK_FLAGS "-openmp")
#TARGET_LINK_LIBRARIES(BarrelContours ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
#SET_TARGET_PROPERTIES(BarrelContours PROPERTIES LINK_FLAGS "-openmp")
# TARGET_LINK_LIBRARIES(LandmarkCentroids ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
# SET_TARGET_PROPERTIES(LandmarkCentroids PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(AxisVectorField ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(AxisVectorField PROPERTIES LINK_FLAGS "-openmp")