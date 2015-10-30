###############################################################################
#
# CMakeLists.txt file for the GridRipper Runtime library
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2013/10/29
#
###############################################################################

# The minimum version of 'cmake' necessary to build this project
cmake_minimum_required (VERSION 2.8.11)
# http://www.cmake.org/cmake/help/v2.8.11/cmake.html#command:cmake_minimum_required


# The project name and version, release numbers
project (Gripper)
set     (Gripper_VERSION_MAJOR 0)
set     (Gripper_VERSION_MINOR 3)
set     (Gripper_VERSION_PATCH 0)
# http://www.cmake.org/cmake/help/v2.8.11/cmake.html#command:project
# http://www.cmake.org/cmake/help/v2.8.11/cmake.html#command:set


# Append our project's include directory to the "#include <>" paths
include_directories ("${PROJECT_SOURCE_DIR}/inc/")
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})
# http://www.cmake.org/cmake/help/v2.8.11/cmake.html#command:include_directories


# Header definitions based on configuration
configure_file (
  "${PROJECT_SOURCE_DIR}/inc/Gripper/Gripper_Config.hpp.in"
  "${PROJECT_SOURCE_DIR}/inc/Gripper/Gripper_Config.hpp"
  )
# http://www.cmake.org/cmake/help/v2.8.11/cmake.html#command:configure_file


# Adding source code files according to configuration
set (Gripper_BUILD) # Variable to hold ALL files to build and be visible in IDE
set (Gripper_HDRS)
set (Gripper_SRCS)

set (Gripper_HDRS ${Gripper_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/Gripper_Config.hpp)
set (Gripper_HDRS ${Gripper_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/Gripper_Export.hpp)

set (Gripper_HDRS ${Gripper_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/Runtime.hpp)
set (Gripper_HDRS ${Gripper_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/Logger.hpp)
set (Gripper_HDRS ${Gripper_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/AnyOption.hpp)

set (Gripper_HDRS ${Gripper_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/MultipoleDefs.hpp)
set (Gripper_HDRS ${Gripper_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/MultipoleTypes.hpp)

set (Gripper_HDRS ${Gripper_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/PDE.hpp)

set (Gripper_SRCS ${Gripper_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/Runtime.cpp)
set (Gripper_SRCS ${Gripper_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/Logger.cpp)
set (Gripper_SRCS ${Gripper_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/AnyOption.cpp)

set (Gripper_SRCS ${Gripper_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/Radial.cpp)
set (Gripper_SRCS ${Gripper_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/Spherical.cpp)
set (Gripper_SRCS ${Gripper_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/Gaunt.cpp)

set (Gripper_BUILD ${Gripper_BUILD} ${Gripper_HDRS} ${Gripper_SRCS})

source_group ("Headers" FILES ${Gripper_HDRS})
source_group ("Sources" FILES ${Gripper_SRCS})

if (GR_USE_STL)

	set (Gripper_STL_HDRS)
	set (Gripper_STL_SRCS)
	
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlRuntime.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlLogger.hpp)
	
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlGnuPlotter.hpp)
	
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlMultipoleDefs.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlMultipoleTypes.hpp)
	
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlGauntIndex.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlGauntCoefficient.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlGauntMatrix.hpp)

	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlRadialIndex.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlRadialExtent.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlRadialVector.hpp)
	
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlSphericalIndex.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlSphericalExtent.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlSphericalVector.hpp)
	
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlExpansionIndex.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlExpansionExtent.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlExpansionConstantFactor.hpp)
	set (Gripper_STL_HDRS ${Gripper_STL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/stl/stlExpansionField.hpp)
	
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlRuntime.cpp)
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlLogger.cpp)
	
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlGauntIndex.cpp)
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlGauntCoefficient.cpp)
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlGauntMatrix.cpp)
	
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlRadialIndex.cpp)
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlRadialExtent.cpp)
	
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlSphericalIndex.cpp)
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlSphericalExtent.cpp)
	
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlExpansionIndex.cpp)
	set (Gripper_STL_SRCS ${Gripper_STL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/stl/stlExpansionExtent.cpp)
	
	set (Gripper_BUILD ${Gripper_BUILD} ${Gripper_STL_HDRS} ${Gripper_STL_SRCS})
	
	source_group ("Headers\\STL" FILES ${Gripper_STL_HDRS})
	source_group ("Sources\\STL" FILES ${Gripper_STL_SRCS})
	
endif (GR_USE_STL)

if (GR_USE_OPENCL)

	set (Gripper_CL_HDRS)
	set (Gripper_CL_SRCS)
	set (Gripper_CL_KERN)
	
	set (Gripper_CL_HDRS ${Gripper_CL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/cl/clRuntime.hpp)
	set (Gripper_CL_HDRS ${Gripper_CL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/cl/clLogger.hpp)
	set (Gripper_CL_HDRS ${Gripper_CL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/cl/clGaunt.hpp)
	
	set (Gripper_CL_HDRS ${Gripper_CL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/cl/clMultipoleDefs.hpp)
	set (Gripper_CL_HDRS ${Gripper_CL_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/cl/clMultipoleTypes.hpp)

	set (Gripper_CL_SRCS ${Gripper_CL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/clRuntime.cpp)
	set (Gripper_CL_SRCS ${Gripper_CL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/clLogger.cpp)
	set (Gripper_CL_SRCS ${Gripper_CL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/clGaunt.cpp)
	
	set (Gripper_CL_SRCS ${Gripper_CL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/clRadial.cpp)
	set (Gripper_CL_SRCS ${Gripper_CL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/clSpherical.cpp)
	set (Gripper_CL_SRCS ${Gripper_CL_SRCS} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/clField.cpp)
	
	set (Gripper_CL_KERN ${Gripper_CL_KERN} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/kernels/Radial.cl)
	set (Gripper_CL_KERN ${Gripper_CL_KERN} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/kernels/Spherical.cl)
	set (Gripper_CL_KERN ${Gripper_CL_KERN} ${PROJECT_SOURCE_DIR}/src/Gripper/cl/kernels/Field.cl)
	
	set (Gripper_BUILD ${Gripper_BUILD} ${Gripper_CL_HDRS} ${Gripper_CL_SRCS} ${Gripper_CL_KERN})
	
	source_group ("Headers\\CL" FILES ${Gripper_CL_HDRS})
	source_group ("Sources\\CL" FILES ${Gripper_CL_SRCS})
	source_group ("Kernels\\CL" FILES ${Gripper_CL_KERN})
	
endif (GR_USE_OPENCL)

if (GR_BUILD_VISUALIZER_SUPPORT)

	set (Gripper_QT_HDRS)
	set (Gripper_QT_SRCS)
	set (Gripper_MOC_SRCS)
	set (Gripper_GLSL)
	
	set (Gripper_QT_HDRS ${Gripper_QT_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/qt/qtRuntime.hpp)
	set (Gripper_QT_HDRS ${Gripper_QT_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/qt/qtLogger.hpp)
	set (Gripper_QT_HDRS ${Gripper_QT_HDRS} ${PROJECT_SOURCE_DIR}/inc/Gripper/qt/qtGaunt.hpp)
	
	set (Gripper_QT_SRCS ${Gripper_QT_SRCS} ${PROJECT_SOURCE_DIR}/inc/Gripper/qt/qtRuntime.cpp)
	set (Gripper_QT_SRCS ${Gripper_QT_SRCS} ${PROJECT_SOURCE_DIR}/inc/Gripper/qt/qtLogger.cpp)
	set (Gripper_QT_SRCS ${Gripper_QT_SRCS} ${PROJECT_SOURCE_DIR}/inc/Gripper/qt/qtGaunt.cpp)
	
	set (Gripper_GLSL ${Gripper_GLSL} ${PROJECT_SOURCE_DIR}/src/Gripper/qt/shaders/Component_view_vertex.glsl)
	set (Gripper_GLSL ${Gripper_GLSL} ${PROJECT_SOURCE_DIR}/src/Gripper/qt/shaders/Component_view_fragment.glsl)
	
	qt_wrap_cpp (Gripper Gripper_MOC_SRCS ${Gripper_QT_HDRS})
	
	set (Gripper_BUILD ${Gripper_BUILD} ${Gripper_QT_HDRS} ${Gripper_QT_SRCS} ${Gripper_MOC_SRCS} ${Gripper_GLSL})
	
	source_group ("Headers\\Qt" FILES ${Gripper_QT_HDRS})
	source_group ("Sources\\Qt" FILES ${Gripper_QT_SRCS})
	source_group ("Sources\\Qt\\Generated\ Files" FILES ${Gripper_MOC_SRCS})
	source_group ("Shaders\\GLSL" FILES ${Gripper_GLSL})
	
endif (GR_BUILD_VISUALIZER_SUPPORT)


# Adding library target for build
add_library (Gripper SHARED ${Gripper_BUILD})
target_link_libraries (Gripper ${Gripper_ADDITIONAL_LIBS})
# http://www.cmake.org/cmake/help/v2.8.11/cmake.html#command:add_library


# Create project groups for IDEs
set_target_properties (Gripper PROPERTIES FOLDER "Library")


# Add platform specific compile options
if(MSVC_VERSION)
#set_target_properties (Gripper PROPERTIES
#                       COMPILE_FLAGS ${MSVC_COMPILER_FLAGS}
#                       LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)