cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(CortexCoordinates)

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

#SET(Sources src/main.cpp src/typedefs.h src/basics.h src/basics.cpp src/amiraReader.h src/amiraReader.cpp src/geometry.h src/geometry.cpp)
SET(Sources ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp src/geometry.h src/geometry.cpp)
SET(Sources2 ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp ../common/barrel_field.h ../common/barrel_field.cpp)

ADD_EXECUTABLE(CortexCoordinates  src/main.cpp ${Sources} )
ADD_EXECUTABLE(SurfaceReconstruction  src/surfaces.cpp ${Sources} )
ADD_EXECUTABLE(ColumnOrientations  src/axis_vectors.cpp ${Sources2} )

#TARGET_LINK_LIBRARIES(CortexCoordinates ITKCommon ITKIO)
# TARGET_LINK_LIBRARIES(CortexCoordinates ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid)
TARGET_LINK_LIBRARIES(CortexCoordinates ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(CortexCoordinates PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(SurfaceReconstruction ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(SurfaceReconstruction PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(ColumnOrientations ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(ColumnOrientations PROPERTIES LINK_FLAGS "-openmp")