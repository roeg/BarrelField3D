cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(NeuroRegistration)
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
# INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(src ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

SET(Sources src/main.cpp src/utilities.h src/utilities.cpp src/morph_reg.h src/morph_reg.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/barrel_field.h ../common/barrel_field.cpp ../common/amiraReader.h ../common/amiraReader.cpp)
SET(Sources2 src/l1_reg.cpp src/utilities.h src/utilities.cpp src/morph_reg.h src/morph_reg.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/barrel_field.h ../common/barrel_field.cpp ../common/amiraReader.h ../common/amiraReader.cpp)
SET(Sources3 src/landmark_reg.cpp src/utilities.h src/utilities.cpp src/morph_reg.h src/morph_reg.cpp ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/barrel_field.h ../common/barrel_field.cpp ../common/amiraReader.h ../common/amiraReader.cpp)

ADD_EXECUTABLE(NeuroRegistration ${Sources} )
ADD_EXECUTABLE(NeuroRegistrationL1 ${Sources2} )
ADD_EXECUTABLE(LandmarkRegistration ${Sources3} )

TARGET_LINK_LIBRARIES(NeuroRegistration ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(NeuroRegistration PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(NeuroRegistrationL1 ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(NeuroRegistrationL1 PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(LandmarkRegistration ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(LandmarkRegistration PROPERTIES LINK_FLAGS "-openmp")