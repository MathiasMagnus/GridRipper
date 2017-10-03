###############################################################################
#
# CMakeLists.txt file for the GridRipper project
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2013/10/29
#
###############################################################################

# The minimum version of 'cmake' necessary to build this project
cmake_minimum_required (VERSION 2.8.11)

# The project name and version, release numbers
project (GridRipper LANGUAGES CXX)
set     (GridRipper_VERSION_MAJOR 0)
set     (GridRipper_VERSION_MINOR 1)
set     (GridRipper_VERSION_PATCH 0)

# Behavioural options for the project
option (GR_USE_STL "Build the STL compute back-end" ON)
option (GR_BUILD_EXAMPLES "Build examples" ON)
option (GR_BUILD_TESTS "Build unit tests" OFF)

#
# Set platform-independent variables
#

  # Append our project's cmake directory to cmake's own module path.
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

  # Use solution folders in IDEs
  set_property (GLOBAL PROPERTY USE_FOLDERS ON)

  # Collect dependencies (let me deal with PackageConfig later)
  find_package(Mathematica REQUIRED)

  # Set configuration-specific variables
  set (Gripper_ADDITIONAL_LIBS)
  set (Gripper_ADDITIONAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/Gripper/inc/)

#
# Set platform-specific variables
#

  # If compiling on Windows
  if(WIN32)
    # The VS 15 fork of CMake contains "MSVC" in its version string
    if(CMAKE_VERSION MATCHES "MSVC")
      message(STATUS "Visual Studio IDE integration detected.")
      message(STATUS "Disabling CMake's automatic reconfiguration in favor of IDE detection.")
      set(CMAKE_SUPPRESS_REGENERATION ON CACHE BOOL "Disables pre-build rule of checking script timestamps")
    endif(CMAKE_VERSION MATCHES "MSVC")
  endif(WIN32)

#
# Set compiler-specific variables
#

  # If compiling with Visual Studio
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

    set (MSVC_CXX_COMPILER_FLAGS)
	set (MSVC_CXX_COMPILER_FLAGS_DEBUG)
	set (MSVC_CXX_COMPILER_FLAGS_RELEASE)
    set (MSVC_CXX_LINKER_FLAGS)

	# Enable C++ exception handling (floating point exceptions in Debug builds)
    set (MSVC_CXX_COMPILER_FLAGS_DEBUG "${MSVC_CXX_COMPILER_FLAGS} /EHa")
    set (MSVC_CXX_COMPILER_FLAGS_RELEASE "${MSVC_CXX_COMPILER_FLAGS} /EHsc /Qpar")

    # Disable incremental linking and manifest files for annoying command line errors
    #set (MSVC_CXX_LINKER_FLAGS " /INCREMENTAL:NO /MANIFEST:NO /NOLOGO ")
  
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MSVC_CXX_COMPILER_FLAGS}")
	set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${MSVC_CXX_COMPILER_FLAGS_DEBUG}")
	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MSVC_CXX_COMPILER_FLAGS_RELEASE}")
	
  endif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

  # If compiling with GCC
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")

    set (GNU_CXX_COMPILER_FLAGS)
    
    # Enable C++11 language and STL support
    set (GNU_CXX_COMPILER_FLAGS "${GNU_CXX_COMPILER_FLAGS} -std=c++11")
    # Enable multi-threaded application
    set (GNU_CXX_COMPILER_FLAGS "${GNU_CXX_COMPILER_FLAGS} -pthread")
    # Enable C++ exception handler support
    set (GNU_CXX_COMPILER_FLAGS "${GNU_CXX_COMPILER_FLAGS} -fexceptions")
    # Enable intrinsic operations native to the machine
    set (GNU_CXX_COMPILER_FLAGS "${GNU_CXX_COMPILER_FLAGS} -march=native")
    # Link statically to std
    set (GNU_CXX_COMPILER_FLAGS "${GNU_CXX_COMPILER_FLAGS} -static-libstdc++")
    
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GNU_CXX_COMPILER_FLAGS}")
  	
  endif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  
  # If compiling with Clang
  if (CMAKE_CXX_COMPLIER_ID MATCHES "Clang")
  
    set (CLANG_CXX_COMPILER_FLAGS)
    
    # Enable C++11 language and STL support
    set (CLANG_CXX_COMPILER_FLAGS "${CLANG_CXX_COMPILER_FLAGS} -std=c++11")
    # Enable Clang's own STL implementation
    set (CLANG_CXX_COMPILER_FLAGS "${CLANG_CXX_COMPILER_FLAGS} -stdlib=libc++")
    
    set(CMAKE_CXX_FLAGS ${CLANG_CXX_COMPILER_FLAGS})
  
  endif (CMAKE_CXX_COMPLIER_ID MATCHES "Clang")

#
# Relay control to subprojects
#

# Process targets inside project
#add_subdirectory (Gripper)

if (GR_BUILD_EXAMPLES)
  add_subdirectory (examples)
endif (GR_BUILD_EXAMPLES)

if (GR_BUILD_TESTS)
  enable_testing()
  add_subdirectory (tests)
endif (GR_BUILD_TESTS)