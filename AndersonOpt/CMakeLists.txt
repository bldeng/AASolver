cmake_minimum_required(VERSION 3.1)

project(AndersonOpt)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

add_executable(ARAP2D
	Types.h
	Constraint.h
	SPDSolver.h
	MeshTypes.h
	SparseMatrixStorage.h
	AndersonAcceleration.h
	GeometrySolver.h
	FileParser.h
    Parameters.h
	ARAP2D.cpp
	)
	
add_executable(ARAP3D
	Types.h
	Constraint.h
	SPDSolver.h
	MeshTypes.h
	SparseMatrixStorage.h
	AndersonAcceleration.h
	GeometrySolver.h
	FileParser.h
	TetMeshIO.h
    Parameters.h
	ARAP3D.cpp
	)	
	
	
add_executable(PlanarityOpt
	Types.h
	Constraint.h
	SPDSolver.h
	MeshTypes.h
	SparseMatrixStorage.h
	AndersonAcceleration.h
	GeometrySolver.h
    Parameters.h
	PlanarityOpt.cpp
	)
	
add_executable(WireMeshOpt
	Types.h
	Constraint.h
	SPDSolver.h
	MeshTypes.h
	SparseMatrixStorage.h
	AndersonAcceleration.h
	GeometrySolver.h
    Parameters.h
	WireMeshOpt.cpp
	)

add_executable(TetMeshSimulation
            Parameters.h
            Types.h
            Constraint.h
            SPDSolver.h
            MeshTypes.h
            SparseMatrixStorage.h
            AndersonAcceleration.h
            DynamicsSolver.h
            SimulationOpt.h
            FileParser.h
            TetMeshIO.h
            TetMeshSimulation.h
            TetMeshSimulation.cpp
    )

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message("Current System: MacOS")
else()
        target_compile_definitions(TetMeshSimulation PUBLIC USE_GL_FOLDER)
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	message("Clang compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	set(OPENMP ON CACHE BOOL "OpenMP")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
	message("AppleClang compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	set(OPENMP OFF CACHE BOOL "OpenMP")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message("GNU compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	set(OPENMP ON CACHE BOOL "OpenMP")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	message("MSVC compiler found.")
	add_definitions(/DUSE_MSVC)
	add_definitions(/D_USE_MATH_DEFINES)
	add_definitions(/DNOMINMAX)
	set(OPENMP ON CACHE BOOL "OpenMP")
endif()




# Add the OpenMesh folder into include path
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/external")

# OpenMesh library
set(OPENMESH_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/OpenMesh/Core")
file(GLOB_RECURSE OPENMESH_HEADERS "${OPENMESH_SOURCE_DIR}/*.hh" "${OPENMESH_SOURCE_DIR}/*.h")
file(GLOB_RECURSE OPENMESH_SOURCES "${OPENMESH_SOURCE_DIR}/*.cc")
add_library(OpenMesh STATIC ${OPENMESH_HEADERS} ${OPENMESH_SOURCES})
add_definitions(-DOM_STATIC_BUILD)

target_link_libraries(ARAP2D OpenMesh)
target_link_libraries(ARAP3D OpenMesh)
target_link_libraries(PlanarityOpt OpenMesh)
target_link_libraries(WireMeshOpt OpenMesh)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/Eigen/Dense)
	message("Found Eigen in local folder")
	set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen")
	target_include_directories(ARAP2D SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(ARAP3D SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(PlanarityOpt SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(WireMeshOpt SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
    target_include_directories(TetMeshSimulation SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
else()
	find_package(Eigen3 REQUIRED)
	if(EIGEN3_FOUND)
		target_include_directories(ARAP2D SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
		target_include_directories(ARAP3D SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
		target_include_directories(PlanarityOpt SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
		target_include_directories(WireMeshOpt SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
        target_include_directories(TetMeshSimulation SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	endif()
endif()


#OpenMP
set(OPENMP ON CACHE BOOL "OpenMP")
if(OPENMP)
  FIND_PACKAGE(OpenMP QUIET)
  if(OPENMP_FOUND)    
      target_compile_options(ARAP2D PUBLIC "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_definitions(ARAP2D PUBLIC "$<$<CONFIG:RELEASE>:USE_OPENMP>")
      target_link_libraries(ARAP2D "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")   
      target_compile_options(ARAP3D PUBLIC "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_definitions(ARAP3D PUBLIC "$<$<CONFIG:RELEASE>:USE_OPENMP>")
      target_link_libraries(ARAP3D "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")   
      target_compile_options(PlanarityOpt PUBLIC "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_definitions(PlanarityOpt PUBLIC "$<$<CONFIG:RELEASE>:USE_OPENMP>")
      target_link_libraries(PlanarityOpt "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")  
      target_compile_options(WireMeshOpt PUBLIC "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_definitions(WireMeshOpt PUBLIC "$<$<CONFIG:RELEASE>:USE_OPENMP>")
      target_link_libraries(WireMeshOpt "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_options(TetMeshSimulation PUBLIC "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_definitions(TetMeshSimulation PUBLIC "$<$<CONFIG:RELEASE>:USE_OPENMP>")
      target_link_libraries(TetMeshSimulation "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      add_definitions(-DEIGEN_DONT_PARALLELIZE)
  else()
      message("OpenMP not found.")
  endif()
endif()

#OpenGL
find_package(OpenGL REQUIRED)
target_include_directories(TetMeshSimulation SYSTEM PUBLIC ${OPENGL_INCLUDE_DIR})

find_package(GLUT REQUIRED)
target_include_directories(TetMeshSimulation SYSTEM PUBLIC ${GLUT_INCLUDE_DIR})

target_link_libraries(TetMeshSimulation OpenMesh ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

if(UNIX AND NOT APPLE)
find_package(Threads REQUIRED)
target_link_libraries(WireMeshOpt ${CMAKE_THREAD_LIBS_INIT})
endif()
