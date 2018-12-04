cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(Connectome)

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

# SET(MapperSources src/singleneuronmain.cpp ../common/spatialgraphset.h ../common/spatialgraphset.cpp)
SET(MatrixVisSources src/connectionMatrixVisualization.cpp)
SET(InnervationHistSources src/innervationHistogram.cpp)
SET(TripletMotifHistSources src/tripletMotifHistogram.cpp)
SET(UniquifyMatrixRowsSources src/uniqueConnectionMatrixRows.cpp)
SET(SynapsesPerCellSources src/synapsesPerCellAnalysis.cpp)
SET(SynapsesAnalysisSources src/synapsesAnalysis.cpp)
SET(GenerateTruncatedCellSources src/generateTruncatedCell.cpp)
SET(HelperSources ../Interneuron/src/helper.h ../Interneuron/src/helper.cpp)
SET(SimpleAxonMoverSources src/moveAxonsSimple.cpp)
SET(TestRealizationSources src/generateTestRealization.cpp)
SET(MatrixSortSources src/sortConnectionMatrix.cpp)
SET(PrePostIDAssignmentSources src/assignPrePostIDs.cpp)
SET(CommonSources ../common/typedefs.h ../common/basics.h ../common/basics.cpp ../common/amiraReader.h ../common/amiraReader.cpp ../common/profile.h ../common/profile.cpp src/matrixanalyzer.h src/matrixanalyzer.cpp ../common/barrel_field.h ../common/barrel_field.cpp)

# ADD_EXECUTABLE(SingleNeuronSynapseMapper ${MapperSources} ${CommonSources})
ADD_EXECUTABLE(MatrixVisualization ${MatrixVisSources} ${CommonSources})
ADD_EXECUTABLE(InnervationHistogram ${InnervationHistSources} ${CommonSources})
ADD_EXECUTABLE(TripletMotifHistogram ${TripletMotifHistSources} ${CommonSources})
ADD_EXECUTABLE(UniquifyMatrixRows ${UniquifyMatrixRowsSources} ${CommonSources})
ADD_EXECUTABLE(SynapsesPerCellAnalysis ${SynapsesPerCellSources} ${CommonSources})
ADD_EXECUTABLE(SynapsesAnalysis ${SynapsesAnalysisSources} ${CommonSources} ${HelperSources})
ADD_EXECUTABLE(GenerateTruncatedCell ${GenerateTruncatedCellSources} ${CommonSources} ${HelperSources})
ADD_EXECUTABLE(SimpleAxonMover ${SimpleAxonMoverSources} ${CommonSources})
ADD_EXECUTABLE(GenerateTestRealization ${TestRealizationSources} ${CommonSources})
ADD_EXECUTABLE(SortConnectionMatrix ${MatrixSortSources} ${CommonSources})
ADD_EXECUTABLE(PrePostIDAssignment ${PrePostIDAssignmentSources} ${CommonSources})

# TARGET_LINK_LIBRARIES(SingleNeuronSynapseMapper ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
# SET_TARGET_PROPERTIES(SingleNeuronSynapseMapper PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(MatrixVisualization ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(MatrixVisualization PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(InnervationHistogram ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(InnervationHistogram PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(TripletMotifHistogram ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(TripletMotifHistogram PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(UniquifyMatrixRows ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(UniquifyMatrixRows PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(SynapsesPerCellAnalysis ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(SynapsesPerCellAnalysis PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(SynapsesAnalysis ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(SynapsesAnalysis PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(GenerateTruncatedCell ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(GenerateTruncatedCell PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(SimpleAxonMover ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(SimpleAxonMover PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(GenerateTestRealization ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(GenerateTestRealization PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(SortConnectionMatrix ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(SortConnectionMatrix PROPERTIES LINK_FLAGS "-openmp")
TARGET_LINK_LIBRARIES(PrePostIDAssignment ITKBasicFilters ITKCommon ITKIO ITKAlgorithms vtkImaging vtkFiltering vtkCommon vtkGraphics vtkIO vtkHybrid ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
SET_TARGET_PROPERTIES(PrePostIDAssignment PROPERTIES LINK_FLAGS "-openmp")